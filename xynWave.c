/*
 * xynWave.c
 *
 *  Created on: Feb 5, 2011
 *      Author: mousomer
 *
 *
 *      symmetry types for symmetric biorthogonal transforms:
 *      half-symmetric: 32100123
 *      whole-symmetric: 3210123
 *       data				lowpass			highpass
 *		length	 head		start	end		start	end
 *		even	lowpass		whole	half	half	whole
 *		even	highpass	half	whole	whole	half
 *		odd		lowpass		whole	whole	half	half
 *		odd		highpass	half	half	whole	whole
 *
 *      via lifting:
 *      	whole symmetric means a_0+=c*(2*a_1)
 *      	half  symmetric means a_0+=c*(a_1+a_0)
 *
 *
 *      algorithmic convection:
 *      			inPlane											non-inPlace
 *      ---------------------------------------------------------------------------------------
 *      lift		inp+=2, out+=2, transform: a 1 a				not implemented
 *
 *      normal		inp++, out[2]+=2, change out-index (oddity)		inp++, out[2]+=1, change oddity
 *
 *
 *
 *      ToDo list:
 *      Align highpass / lowpass filter positions using filt->pos during fwd/inv filter transform.
 *
 */


#include "libxyn.h"
#define LOW_LEN(len, oddity) len/2 + (len%2 * (1-oddity))

static int FixWaveShift = 10;
static int FixWaveScale = 1000;


/** declarations of static functions **/
/* filter generators (filter are global variables) */
static void wavInitHaar();
static void wavInitCDF53();
static void wavInitCDF97();
static void wavInitTest();

static void wavInitLiftCDF53();
static void wavInitLiftCDF97();

/*static void wavInitPloyHaar();
static void wavInitPolyCDF53();
static void wavInitPolyCDF97();*/

static void genPolyFilterFromWave ( waveFilterPtr wfp, wavePolyPtr pfp );
/** QMF symmetries and filter operations **/
static void wavCompleteWavFilterData ( waveFilterPtr wfp );
/** QMF auxiliaries: **/
static int wavFilterQMF_f	( float* out, float* inp, int len, int zeroPos );
static int wavFilterFold_f	( float (*out)[FILTER_LEN], int* outLen, float* inp, int len, int zeroPos );
static int wavFilterQMF_i	( int* out, int* inp, int len, int zeroPos );
static int wavFilterSym_f	( float* data, int len, int zeroPos );
static int wavCoeffsIntFromFloat ( int* int_coeffs, float* float_coeffs, int n );

/* lift stage in-place */
static int wavLiftingStageInPlace_f ( float* data, long len, int start, register float a );
static int wavLiftingStageInPlace_i (   int* data, long len, int start, register int a );
static int wavLiftingStageInPlace_d ( double* data, long len, int start, register double a );
static int wavLiftingStageInPlace_fVect ( float* data, long vLen, int vSize, int start, register float a );
static int wavLiftingStageInPlace_iVect (   int* data, long vLen, int vSize, int start, register int a );
static int wavLiftingStageInPlace_dVect ( double* data, long vLen, int vSize, int start, register double a );
/* filter transform - main loops */
static int wavTransformHiLow_f		( float* dataInp[2], int ji, float* dataOut[2], int jo, long len, bool oddity, waveFilterPtr filt );
static int wavInvTransformHiLow_f	( float* dataInp[2], int ji, float* dataOut, int jo, long len, bool oddity, waveFilterPtr filt );

/* lift stage not in-place
static int wavLiftingStage_f ( float* inpData, float* outData, long len, int start, register float a );
static int wavLiftingStage_i (   int* inpData,   int* outData, long len, int start, register int a );
static int wavLiftingStage_fVect ( float* inpData, float* outData, long vLen, int vSize, int start, register float a );
static int wavLiftingStage_iVect (   int* inpData,   int* outData, long vLen, int vSize, int start, register int a );
*/

/* packing / unpacking (reordering after/before lifting) */
static int wavPack1d_f		( float* data,  long len, int oddity );
static int wavPack1d_i		(   int* data,  long len, int oddity );
static int wavPack1d_d 		( double* data, long len, int oddity );
static int wavUnpack1d_f 	( float* data,  long len, int oddity );
static int wavUnpack1d_i 	(   int* data,  long len, int oddity );

/* calculate last lowpass position */
#define SCALE_END(len, oddity) ( (len%2)&&oddity ) || ( (1-len%2)&&(1-oddity) )

/*   global wavelet transform filters */
struct waveletFilter wavHaar, wavCDF97, wavCDF53, wavTest;

struct wavePolyFilter wavPolyHaar, wavPolyCDF97, wavPolyCDF53;

struct waveLiftFilter wavLiftHaar, wavLiftCDF97, wavLiftCDF53;

/** filter print **/
int wavPrint ( waveFilterPtr filt, FILE* fp )
{
	fprintf (fp, "wavelet filter %s\n", filt->name);

	fprintf( fp, "lowPass analysis   (start at %d): ", filt->pos_anal[0]);
	displayArray_f ( filt->f_anal[0], filt->len_anal[0], fp );

	fprintf( fp, "highPass analysis  (start at %d): ", filt->pos_anal[1]);
	displayArray_f  ( filt->f_anal[1], filt->len_anal[1], fp );

	fprintf( fp, "lowPass synthesis  (start at %d): ", filt->pos_synt[0]);
	displayArray_f  ( filt->f_synt[0], filt->len_synt[0], fp );

	fprintf( fp, "highPass synthesis (start at %d): ", filt->pos_synt[1]);
	displayArray_f  ( filt->f_synt[1], filt->len_synt[1], fp );

	fprintf( fp, "\nfolded elements: lowpass analysis\t");
	displayArray_i  ( filt->len_anal_folded[0], filt->len_anal[0], fp ); NL;
	for (int i=0; i<filt->len_anal[0]; i++)
		displayArray_f  ( filt->f_anal_folded[0][i], filt->len_anal_folded[0][i], fp );

	fprintf( fp, "\nfolded elements: highpass analysis\t");
	displayArray_i  ( filt->len_anal_folded[1], filt->len_anal[1], fp ); NL;
	for (int i=0; i<filt->len_anal[1]; i++)
		displayArray_f  ( filt->f_anal_folded[1][i], filt->len_anal_folded[1][i], fp );

	fprintf( fp, "\nfolded elements: lowpass synthesis\t");
	displayArray_i  ( filt->len_synt_folded[0], filt->len_synt[0], fp ); NL;
	for (int i=0; i<filt->len_synt[0]; i++)
		displayArray_f  ( filt->f_synt_folded[0][i], filt->len_synt_folded[0][i], fp );

	fprintf( fp, "\nfolded elements: highpass synthesis\t");
	displayArray_i  ( filt->len_synt_folded[1], filt->len_synt[1], fp ); NL;
	for (int i=0; i<filt->len_synt[1]; i++)
		displayArray_f  ( filt->f_synt_folded[1][i], filt->len_synt_folded[1][i], fp );
	NL;

	fprintf( fp, "\nfolded elements: lowpass analysis\t");
	displayArray_i  ( filt->len_anal_folded[0], filt->len_anal[0], fp ); NL;
	for (int i=0; i<filt->len_anal[0]; i++)
		displayArray_i  ( filt->i_anal_folded[0][i], filt->len_anal_folded[0][i], fp );

	fprintf( fp, "\nfolded elements: highpass analysis\t");
	displayArray_i  ( filt->len_anal_folded[1], filt->len_anal[1], fp ); NL;
	for (int i=0; i<filt->len_anal[1]; i++)
		displayArray_i  ( filt->i_anal_folded[1][i], filt->len_anal_folded[1][i], fp );

	fprintf( fp, "\nfolded elements: lowpass synthesis\t");
	displayArray_i  ( filt->len_synt_folded[0], filt->len_synt[0], fp ); NL;
	for (int i=0; i<filt->len_synt[0]; i++)
		displayArray_i  ( filt->i_synt_folded[0][i], filt->len_synt_folded[0][i], fp );

	fprintf( fp, "\nfolded elements: highpass synthesis\t");
	displayArray_i  ( filt->len_synt_folded[1], filt->len_synt[1], fp ); NL;
	for (int i=0; i<filt->len_synt[1]; i++)
		displayArray_i  ( filt->i_synt_folded[1][i], filt->len_synt_folded[1][i], fp );
	NL;
	return OK;
}

int wavPolyPrint ( wavePolyPtr filt, FILE* fp )
{
	fprintf (fp, "wavelet filter %s\n", filt->name);

	fprintf( fp, "lowPass analysis   (start at %d/%d): ", filt->pos_anal[0][0], filt->pos_anal[0][1]);
	displayArray_f ( filt->f_anal[0][0], filt->len_anal[0][0], fp );
	displayArray_f ( filt->f_anal[0][1], filt->len_anal[0][1], fp );

	fprintf( fp, "highPass analysis  (start at %d/%d): ", filt->pos_anal[1][0], filt->pos_anal[1][1]);
	displayArray_f  ( filt->f_anal[1][0], filt->len_anal[1][0], fp );
	displayArray_f  ( filt->f_anal[1][1], filt->len_anal[1][1], fp );

	fprintf( fp, "lowPass synthesis  (start at %d/%d): ", filt->pos_synt[0][0], filt->pos_synt[0][1]);
	displayArray_f  ( filt->f_synt[0][0], filt->len_synt[0][0], fp );
	displayArray_f  ( filt->f_synt[0][1], filt->len_synt[0][1], fp );

	fprintf( fp, "highPass synthesis (start at %d/%d): ", filt->pos_synt[1][0], filt->pos_synt[1][1]);
	displayArray_f  ( filt->f_synt[1][0], filt->len_synt[1][0], fp );
	displayArray_f  ( filt->f_synt[1][1], filt->len_synt[1][1], fp );


	fprintf( fp, "\nfolded elements:\n");

	for (int p=0; p<=1; p++)
		for (int o=0; o<=1; o++)
		{
			fprintf( fp, "\n analysis %d/%d lengths ", p, o);
			displayArray_i  ( filt->len_anal_folded[p][o], filt->len_anal[p][o], fp ); NL;
			for (int i=0; i<filt->len_anal[p][o]; i++)
				displayArray_f  ( filt->f_anal_folded[p][o][i], filt->len_anal_folded[p][o][i], fp );

			fprintf( fp, "\n synthesis %d/%d lengths ", p, o);
			displayArray_i  ( filt->len_synt_folded[p][o], filt->len_synt[p][o], fp ); NL;
			for (int i=0; i<filt->len_synt[p][o]; i++)
				displayArray_f  ( filt->f_synt_folded[p][o][i], filt->len_synt_folded[p][o][i], fp );
		}

	for (int p=0; p<=1; p++)
		for (int o=0; o<=1; o++)
		{
			fprintf( fp, "\n analysis %d/%d lengths ", p, o);
			displayArray_i  ( filt->len_anal_folded[p][o], filt->len_anal[p][o], fp ); NL;
			for (int i=0; i<filt->len_anal[p][o]; i++)
				displayArray_i  ( filt->i_anal_folded[p][o][i], filt->len_anal_folded[p][o][i], fp );

			fprintf( fp, "\n synthesis %d/%d lengths ", p, o);
			displayArray_i  ( filt->len_synt_folded[p][o], filt->len_synt[p][o], fp ); NL;
			for (int i=0; i<filt->len_synt[p][o]; i++)
				displayArray_i  ( filt->i_synt_folded[p][o][i], filt->len_synt_folded[p][o][i], fp );
		}
	return OK;
}


int wavLiftPrint ( waveLiftPtr filt, FILE* fp )
{
	fprintf (fp, "wavelet lifting filter %s, number of predict/update stages: %d\n", filt->name, filt->nPredict);
	fprintf (fp, "wavelet lifting filter %s\n", filt->name);

	fprintf( fp, "floating point update filter: " );
	displayArray_f ( filt->coeff_f[0], filt->nPredict, fp );

	fprintf( fp, "floating point predict filter: " );
	displayArray_f ( filt->coeff_f[1], filt->nPredict, fp );

	fprintf( fp, "floating point lowpass  scaler: %f\n", filt->scaler_f[0] );
	fprintf( fp, "floating point highpass scaler: %f\n", filt->scaler_f[1] );
	fprintf( fp, "integer lowpass  scaler: %d\n", filt->scaler_i[0] );
	fprintf( fp, "integer highpass scaler: %d\n", filt->scaler_i[1] );

	fprintf( fp, "integer update filter: " );
	displayArray_i ( filt->coeff_i[0], filt->nPredict, fp );

	fprintf( fp, "integer predict filter: " );
	displayArray_i ( filt->coeff_i[1], filt->nPredict, fp );


	return OK;
}

/** filter generation by wavelet name **/
waveFilterPtr wavGetWaveletName( const char* name )
{
	waveFilterPtr wp=NULL;
	//printf ("function %s. %s\n", __FUNCTION__, name );
	if ( strncmp(name, "CDF-5-3", 7) == 0 )
	{
		wp = &wavCDF53;
		wavInitCDF53();
		//printf("wavelet: CDF-5-3\n");
	}
	else if ( strncmp(name,"Haar", 4) == 0 )
	{
		wp = &wavHaar;
		wavInitHaar();
		//printf("wavelet: Haar\n");
	}
	else if ( strncmp(name, "CDF-9-7", 7) ==0 )
	{
		wp = &wavCDF97;
		wavInitCDF97();
		//printf("wavelet: CDF-9-7\n");
	}
	else if ( strncmp(name, "test", 4) ==0 )
	{
		wp = &wavTest;
		wavInitTest();
		//printf("wavelet: test\n");
	}
	else
		perror ("no such wavelet filter implemented here.\n");

	return wp;
}


waveLiftPtr wavGetLiftWaveName( const char* name )
{
	waveLiftPtr wp = NULL;
	//printf ("function %s. %s\n", __FUNCTION__, name );
	if ( strncmp(name, "CDF-5-3", 7) == 0 )
	{
		wp = &wavLiftCDF53;
		wavInitLiftCDF53();
	}
	else if ( strncmp(name,"Haar", 4) == 0 )
	{
		ErrPrintMsg("not a valid lifting filter/n");
		//wp = &wavLiftHaar;
		//wavInitHaar();
	}
	else if ( strncmp(name, "CDF-9-7", 7) ==0 )
	{
		wp = &wavLiftCDF97;
		wavInitLiftCDF97();
	}
	return wp;
}

wavePolyPtr wavGetPolyWaveName( const char* name )
{
	waveFilterPtr fp = NULL;
	wavePolyPtr pfp = NULL;
	//printf ("function %s. %s\n", __FUNCTION__, name );
	if ( strncmp(name, "CDF-5-3", 7) == 0 )
	{
		fp = &wavCDF53;
		pfp = &wavPolyCDF53;
		wavInitCDF53();
	}
	else if ( strncmp(name,"Haar", 4) == 0 )
	{
		fp = &wavHaar;
		pfp = &wavPolyHaar;
		wavInitHaar();
	}
	else if ( strncmp(name, "CDF-9-7", 7) ==0 )
	{
		fp = &wavCDF97;
		pfp = &wavPolyCDF97;
		wavInitCDF97();
	}

	genPolyFilterFromWave ( fp, pfp );
	return pfp;
}


dwtInfPtr wavGenDWTInfo ( const char* wavName, enum waveType wtype, bool oddity, bool inPlace )
{
	assert ( wavName );
	dwtInfPtr ip = (dwtInfPtr) malloc (sizeof(struct dwtInfo));
	assert ( ip );

	ip->oddity = (oddity>0) ? 1 : 0;
	ip->inPlace = (inPlace>0) ? 1 : 0;
	ip->wType = wtype;
	ip->fwp = NULL;
	ip->lwp = NULL;

	if ( ip->wType == xynWaveleType )
		ip->fwp = wavGetWaveletName ( wavName );
	else if ( ip->wType == xynLiftingType )
		ip->lwp = wavGetLiftWaveName ( wavName );
	else
		perror ("wavelet type not implemented yet.\n");

	return ip;
}



/**** DWT transforms types:
 * wavTransform 	- normal filter
 * invTransform 	- inverse
 * liftTransform 	- lifting implementation
 * inPlace 			- leave data in place (data -> low, high, low, high.... )
 * HiLow 			- fwd transform divided into two output (fwd) / input (inv) arrays
 * ..._i			- integer transform
 * ..._f			- floating point transform
 * IO:
 * int/float *dataInp (sometimes *inpLow, *inpHigh) 	- input array[s]
 * int/float *dataOut (sometimes *outLow, *outHigh) 	- output array[s]
 * long len (sometimes lenLow, lenHigh)					- array[s] length
 * boolean oddity 		- 0 start with lowpass, 1 start with highpass.
 * boolean inPlave		- tells lifting implementation whether to pack/unpack data
 * boolean rescale		- tells lifting implementation whether to rescale (0) or not (1) the data
 * filt					- DWT filter (waveFilter for normal transform, waveLift for lifting)
 * ****/


/* floating point transform - get 1 output data set, divide it into high half and low half */
int wavTransform_f( float* dataInp, float* dataOut, long len, bool oddity, bool inPlace, waveFilterPtr filt )
{
	assert(dataInp!=NULL); assert(dataOut!=NULL); assert(filt!=NULL);

	float *cpyOut[2] = { dataOut, dataOut };			// work with IO pointer copies
	float *cpyInp[2] = { dataInp, dataInp };

	int jumpOut = 2, jumpInp = 2;
	if (inPlace)
	{
		cpyOut[0] += oddity;
		cpyOut[1] += 1-oddity;
	}
	else
	{
		cpyOut[1] += LOW_LEN(len, oddity);
		jumpOut = 1;
	}

	XYN_DPRINTF(DEBUG_HEADER, "%s %s: len %d; oddity %d. jumps: in %d out %d\n", __FUNCTION__, filt->name, len, oddity, jumpInp, jumpOut);

	wavTransformHiLow_f ( cpyInp, jumpInp, cpyOut, jumpOut, len, oddity, filt );

	return OK;
}

int wavInvTransform_f( float* dataInp, float* dataOut, long len, bool oddity, bool inPlace, waveFilterPtr filt )
{
	assert(dataInp!=NULL); assert(dataOut!=NULL); assert(filt!=NULL);

	float *cpyInp[2] = { dataInp, dataInp };

	int jumpOut = 1, jumpInp = 2;
	if (inPlace)
	{
		cpyInp[0] += oddity;
		cpyInp[1] += 1-oddity;
	}
	else
	{
		cpyInp[1] += LOW_LEN(len, oddity);
		jumpInp = 1;
		jumpOut = 1;
	}

	XYN_DPRINTF(DEBUG_HEADER, "%s %s: len %d; oddity %d. jumps: in %d out %d\t lowLen %ld (%f)\n",
			__FUNCTION__, filt->name, len, oddity, jumpInp, jumpOut, len/2 + (len%2 * (1-oddity)), *cpyInp[1] );

	wavInvTransformHiLow_f ( cpyInp, jumpInp, dataOut, jumpOut, len, oddity, filt );

	return OK;
}

/* floating point transform - normal filter, input (1) -> output (2 arrays)
 * set starting positions, boundary lengths, etc.
 * run transform on lower boundary - use folded filters
 * (set fold position)
 * run transform normally - keep changing oddity (from start-> high to low to high...)
 * output[low/high] = filter[low/high] * data
 * run transform on upper
 *
 * Assumptions:
 * 			- filter length (and hence, folded boundaries of the filter) should be shorter than any data buffer
 * 			- filter is symmetric (in form, if not in content) around filt->zeroPos
 *
 * jumpInp should be 2
 * jumpOut  is 2 for inPlace (l h l h...)
 * 			is 1 for packed transform (lllll... hhhhhh...)
 * */
static int wavTransformHiLow_f( float* dataInp[2], int ji, float* dataOut[2], int jo,
							long len, bool oddity, waveFilterPtr filt )
{
	/* check input integrity */
	assert(dataInp[0]!=NULL); assert(dataInp[1]!=NULL); assert(dataOut[0]!=NULL); assert(dataOut[1]!=NULL); assert(filt!=NULL);
	if ( len<XYN_MAX(filt->len_anal[0], filt->len_anal[1] ) )
		ErrMsg("data no long enough for filter\n");

	int fold[2]={0,0}; fold[1-oddity]=1;		// fold[0/1] is lowpass/highpass folded filter index

	XYN_DPRINTF(DEBUG_HEADER,"%s: (len %d). oddity: %d. Folded filters %d/%d\n", __FUNCTION__, len, oddity, fold[0], fold[1]);

	for (long i=0; i<len; i++, oddity = 1-oddity)		// i from 1 - this is second pass type
	{
		if (i < filt->pos_anal[oddity])					// inside lower boundary - call folded filter
		{
			*(dataOut[oddity]) = dotProd_f ( filt->f_anal_folded[oddity][fold[oddity]], dataInp[oddity],
					filt->len_anal_folded[oddity][fold[oddity]] );

			XYN_DPRINTF(DEBUG_DETAIL, "i=%6ld; fold %4d=%7.2f(%d)\t", i, fold[oddity], *dataOut[oddity], oddity);
			//displayArray_f(filt->f_anal_folded[oddity][fold[oddity]], filt->len_anal_folded[oddity][fold[oddity]], stdout);
			fold[oddity]+=2;

			if ( fold[oddity] > filt->pos_anal[oddity] )
				dataInp[oddity] ++;
		}
		else if ( i >= len-filt->pos_anal[oddity] )			// inside upper boundary - call upper filter
		{
			int cFold = filt->len_anal[oddity] - (len-i);	/* calculate current fold (oddity may cancel dependence on lower fold) */
			*(dataOut[oddity]) = dotProd_f ( filt->f_anal_folded[oddity][cFold], dataInp[oddity],
					filt->len_anal_folded[oddity][cFold] );

			XYN_DPRINTF(DEBUG_DETAIL, "\ni=%6ld; fold %4d=%7.2f(%d)\t", i, cFold, *dataOut[oddity], oddity);
			//displayArray_f(filt->f_anal_folded[oddity][cFold], filt->len_anal_folded[oddity][cFold], stdout);

			dataInp[oddity] += ji;
		}
		else
		{
			*(dataOut[oddity]) = dotProd_f ( filt->f_anal[oddity], dataInp[oddity], filt->len_anal[oddity] );

			XYN_DPRINTF(DEBUG_DETAIL, "%4ld=%5.2f(%d)   ", i, *dataOut[oddity], oddity);
			if (i%10==9) XYN_DPRINTF(DEBUG_DETAIL, "\n");
			dataInp[oddity] += ji;
		}

		dataOut[oddity] += jo;								// output pointer advance
	}

	return OK;
}

/* floating point transform - normal filter, input (1) -> output (2 arrays)
 * set starting positions, boundary lengths, etc.
 * run transform on lower boundary - use folded filters
 * (set fold position)
 * run transform normally - keep changing oddity (from start-> high to low to high...)
 * output[low/high] = filter[low/high] * data
 * run transform on upper
 *
 * Assumptions: filter length (and hence, folded boundaries of the filter) should be shorter than any data buffer
 *
 *	jumpInp (ji) is the input jump:
 *				- for packed transform (llll...hhhh...   -> d d d d d...)
 *				- for inPlace transform (lhlhlhlh... -> d d d d d...)
 *	jumpOut (jo) is the output jump - should be 1.
 *
 * lower boundary depends on oddity (0 oddity for starting lowpass)
 * upper boundary depends on end_oddity
 * 						0 for ending lowpss   {oddity=0 len odd, or oddity=1 len even}
 * 						1 for ending highpass {oddity=1 len odd, or oddity=0 len even}
 * 						i.e. end_oddity =  (oddity &  len%2 ) | (1-oddity)&(1-len%2)
 *	algorithm:
 *			calc lowpass (with movement by oddity)
 *			add highpass
 * */
static int wavInvTransformHiLow_f( float* dataInp[2], int ji, float* dataOut, int jo,
									long len, bool oddity, waveFilterPtr filt )
{
	/* check input integrity */
	assert(dataInp[0]!=NULL); assert(dataInp[1]!=NULL); assert(dataOut!=NULL); assert(filt!=NULL);
	if ( len<XYN_MAX(filt->len_synt[0], filt->len_synt[1] ) )
		ErrMsg("data no long enough for filter\n");

	bool od0=oddity, od1=1-oddity;
//	bool end_oddity = (od0 * len%2 ) | ( od1*(1 - len%2) );
	/* start points, lengths. low/high * even(zero)/odd
	 * e.g. 9-7 inverse is 7-9 (l-h) ->
	 * positions are ordinary +-3/+-4 ->
	 * polyphase 	[-2 0 2 / -3 -1 1 3]/[-4 -2 0 2 4 / -3 -1 1 3]
	 * into 	 	[-1 0 1 / -2 -1 0 1 +od] / [-2 -1 0 1 2 / -2 -1 0 1 +od]
	 * lengths: (3/4) (low e/o) / (5/4) (high e/o)
	 * positions: (1/-2 (+oddity) ) / (2/-1 (-oddity) ) */
	int strt[2][2] = { { filt->pos_synt[0]%2, (filt->pos_synt[0]-1)%2 } , { (filt->pos_synt[1])%2, (filt->pos_synt[1]+1)%2 } };
	int lens[2][2] = { { filt->len_synt[0]/2, (filt->len_synt[0]+1)/2 } , { (filt->len_synt[1]+1)/2, (filt->len_synt[1]+1)/2 } };

	XYN_DPRINTF(DEBUG_HEADER,"%s: (len %d). oddity %d. inp/out Jumps %d/%d\ndata input: %p/%p, output: %p/%p\n",
			__FUNCTION__, len, oddity, ji, jo, dataInp[0], dataInp[1], dataOut[0], dataOut[1] );
	XYN_DPRINTF(DEBUG_HEADER,"lengths (l_even/l_odd %d/%d) (h_even/h_odd %d/%d)\npositions (l_even/l_odd %d/%d) (h_even/h_odd %d/%d)\n",
			lens[0][0], lens[0][1], lens[1][0], lens[1][1], strt[0][0], strt[0][1], strt[1][0], strt[1][1] );


	for (long i=0; i<len; i++, od0=od1, od1=1-od1)		// i from 1 - this is second pass type
	{
		/* start lowpass */
		if (i < filt->pos_synt[0])						// inside lower boundary - call folded filter
		{
			XYN_DPRINTF(DEBUG_DETAIL, "fold %d (%d): \t", i, od0);
			//displayArrayJump_f ( filt->f_synt_folded[0][i]+oddity, 2, (filt->len_synt_folded[0][i]+1)/2, stdout);
			//displayArrayJump_f ( dataInp[0], 1, (filt->len_synt_folded[0][i]+1)/2, stdout);
			*dataOut = dotProdJump_f ( filt->f_synt_folded[0][i]+oddity, 2, dataInp[0], 1, (filt->len_synt_folded[0][i]+1)/2 );
		}
		else if ( i >= len-filt->pos_synt[0] )			// inside upper boundary - call upper filter
		{
			int cFold = filt->len_synt[0] - (len-i);	// calculate current fold (od0 may cancel dependence on lower fold)
			bool end_oddity = ( ((i-filt->pos_synt[0])%2) != oddity );
			float* lowFold = filt->f_synt_folded[0][cFold] ;
			XYN_DPRINTF(DEBUG_DETAIL, "low  fold %d (i=%d, end-oddity %d): \n", cFold, i, end_oddity);
			*dataOut = dotProdJump_f ( lowFold+end_oddity, 2, dataInp[0], 1, (filt->len_synt_folded[0][cFold]+1)/2 );
			if (od0==1) dataInp[0] += ji;
		}
		else
		{
			*dataOut = dotProdJump_f ( filt->f_synt[0]+strt[0][od0], 2, dataInp[0], 1, lens[0][od0] );
			//DPRINTF(DEBUG_DETAIL, "%4ld=%7.2f(%d+%d)", i, *dataOut, od0, strt[0][od0]);
			if (od0==1) dataInp[0] += ji;
		}

		//DPRINTF(DEBUG_DETAIL, "%4ld=%7.2f(%d) ", i, *dataOut[od1], od1);
		/* add highpass */
		if (i < filt->pos_synt[1])					// inside lower boundary - call folded filter
		{
			XYN_DPRINTF(DEBUG_DETAIL, "fold %d (%d): \t", i, od1);
			//displayArrayJump_f ( filt->f_synt_folded[1][i]+1-oddity, 2, (filt->len_synt_folded[1][i]+1)/2, stdout);
			//displayArrayJump_f ( dataInp[1], 1, (filt->len_synt_folded[1][i]+1)/2, stdout);
			*dataOut += dotProdJump_f ( filt->f_synt_folded[1][i]+1-oddity, 2, dataInp[1], 1, (filt->len_synt_folded[1][i]+1)/2 );
		}
		else if ( i >= len-filt->pos_synt[1] )			// inside upper boundary - call upper filter
		{
			int cFold = filt->len_synt[1] - (len-i);	// calculate current fold (od1 may cancel dependence on lower fold)

			bool end_oddity = ( ( (i-filt->pos_synt[1])%2) == oddity );
			float* highFold = filt->f_synt_folded[1][cFold];

			XYN_DPRINTF(DEBUG_DETAIL, "\nhigh fold %d (i=%d, end-oddity %d): \n", cFold, i, 1-end_oddity);
			displayArrayJump_f ( highFold+end_oddity, 2, (filt->len_synt_folded[1][cFold]+1)/2, stdout);
			displayArrayJump_f ( dataInp[1], 1, (filt->len_synt_folded[1][cFold]+1)/2, stdout);
			*dataOut += dotProdJump_f ( highFold+end_oddity, 2, dataInp[1], 1, (filt->len_synt_folded[1][cFold]+1)/2 );
			if (od1==0) dataInp[1] += ji;
		}
		else
		{
			*dataOut += dotProdJump_f ( filt->f_synt[1]+strt[1][od1], 2, dataInp[1], 1, lens[1][od1] );
			//DPRINTF(DEBUG_DETAIL, "+(%d+%d)=%7.2f  ", od1, strt[1][od1], *dataOut); if (i%5==4) DPRINTF(DEBUG_DETAIL,"\n");
			if (od1==0) dataInp[1] += ji;
		}

		dataOut++;							// output pointer advance
	}

	return OK;
}


/* normal DWT filter transforms and inverse transforms - integer versions */
int wavTransform_i( int* dataInp, int* dataOut, long len, bool oddity, waveFilterPtr filt )
{
	assert(dataInp!=NULL); assert(filt!=NULL); assert(dataOut!=NULL); assert(len>=0);
	int start[2] = { filt->pos_anal[0]+oddity, filt->pos_anal[1]+1-oddity };
	long startInp = XYN_MIN(start[0], start[1]);
	long endInp = len - startInp;

	int* cpyInp = dataInp;
	int* cpyOut[2] = { dataOut, dataOut + len/2 + oddity*(len%2) };
	XYN_DPRINTF(DEBUG_HEADER, "%s %s: len %d, startLow %d startHigh %d oddity %d\n",
			__FUNCTION__, filt->name, len, start[0], start[1], oddity);

	cpyOut[0] += start[0];
	cpyOut[1] += start[1];
	for ( int i=startInp; i<endInp; i++)
	{
		*cpyOut[oddity]++ = dotProd_i ( filt->i_anal[oddity], cpyInp++, filt->len_anal[oddity] );
		oddity = 1-oddity;
	}

	return OK;
}


int wavInvTransform_i( int* dataInp, int* dataOut, long len, bool oddity, waveFilterPtr filt )
{
	int start = XYN_MAX( filt->pos_synt[0]+oddity, filt->pos_synt[1]+1-oddity);
	int end = len - start;

	int* cpyInp = dataInp;
	int* cpyOut = dataOut;
	XYN_DPRINTF(DEBUG_HEADER, "%s %s: len %d, start %d oddity %d\n",
			__FUNCTION__, filt->name, len, start, oddity);

	cpyOut += start;
	for ( int i=start; i<end; i++)
	{
		*cpyOut++ = dotProd_i ( filt->i_synt[oddity], cpyInp++, filt->len_synt[oddity] );
		oddity = 1-oddity;
	}

	return OK;
}


/** DWT inplace  **/
/* in place transform */
int wavTransformInPlace_i( int* data, long len, bool oddity, waveFilterPtr filt )
{
	int start[2] = { filt->pos_anal[0]+oddity, filt->pos_anal[1]+1-oddity };
	long startInp = XYN_MIN(start[0], start[1]);
	long endInp = len - startInp;

	int* cpyInp = data;
	int* cpyOut[2] = { data, data+ len/2 + oddity*(len%2) };

	cpyOut[0] += start[0];
	cpyOut[1] += start[1];
	for ( int i=startInp; i<endInp; i++)
	{
		*cpyOut[oddity]++ = dotProd_i ( filt->i_anal[oddity], cpyInp++, filt->len_anal[oddity] );
		oddity = 1-oddity;
	}


	return OK;
}

int wavInvTransformInPlace_i( int* data, long len, bool oddity, waveFilterPtr filt )
{
	int start = XYN_MAX( filt->pos_synt[0]+oddity, filt->pos_synt[1]+1-oddity);
	int end = len - start;

	int* cpyInp = data;
	int* cpyOut = data;

	cpyOut += start;
	for ( int i=start; i<end; i++, oddity = 1-oddity)
		*cpyOut++ = dotProd_i ( filt->i_synt[oddity], cpyInp++, filt->len_synt[oddity] );


	return OK;
}



/*** lifting implementation -
 * 			(predict high from low, update low from high) X num_of_predictions
 * 			scale (optional)
 * 			pack (optional)
 *  */
int liftTransform_d( double* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(len, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) (end: %d %d)\n",	__FUNCTION__, data, len, oddity, filt->name, sclEnd, wavEnd);

	for (int k=0; k<filt->nPredict; k++)
	{
		// predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_d ( data, len, oddity, filt->coeff_f[0][k] );

		// update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"update: \t");
		wavLiftingStageInPlace_d ( data, len, 1-oddity, filt->coeff_f[1][k] );
	}

	if ( !noRescale )
	{
		XYN_DPRINTF(DEBUG_PROCESS, "%s rescaling low by %f high by %f\n", __FUNCTION__, filt->scaler_f[1], filt->scaler_f[0] );
		mulVectScalJump_d ( data+oddity, filt->scaler_f[0], len-oddity , 2 );
		mulVectScalJump_d ( data+1-oddity, filt->scaler_f[1], len-1+oddity , 2 );
	}

	// packing
	if ( !inPlace )
	{
		XYN_DPRINTF(DEBUG_PROCESS, "%s packing\n", __FUNCTION__ );
		wavPack1d_d ( data, len, oddity );
	}


	return OK;
}

int liftTransform_f( float* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(len, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) (end: %d %d)\n",
			__FUNCTION__, data, len, oddity, filt->name, sclEnd, wavEnd);

	for (int k=0; k<filt->nPredict; k++)
	{
		// predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_f ( data, len, oddity, filt->coeff_f[0][k] );

		// update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"update: \t");
		wavLiftingStageInPlace_f ( data, len, 1-oddity, filt->coeff_f[1][k] );
	}

	if ( !noRescale )
	{
		XYN_DPRINTF(DEBUG_PROCESS, "%s rescaling low by %f high by %f\n", __FUNCTION__, filt->scaler_f[1], filt->scaler_f[0] );
		mulVectScalJump_f ( data+oddity, filt->scaler_f[0], len-oddity , 2 );
		mulVectScalJump_f ( data+1-oddity, filt->scaler_f[1], len-1+oddity , 2 );
	}

	// packing
	if ( !inPlace )
	{
		XYN_DPRINTF(DEBUG_PROCESS, "%s packing\n", __FUNCTION__ );
		wavPack1d_f ( data, len, oddity );
	}


	return OK;
}


int liftInvTransform_f( float* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(len, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) end: (%d %d)\n",
			__FUNCTION__, data, len, oddity, filt->name, sclEnd, wavEnd);

	// unPack
	if ( !inPlace )
	{
		XYN_DPRINTF(DEBUG_PROCESS, "%s unPacking n", __FUNCTION__ );
		wavUnpack1d_f ( data, len, oddity );
	}

	if ( !noRescale )
	{
		XYN_DPRINTF(DEBUG_PROCESS, "%s rescaling low by %f high by %f\n", __FUNCTION__, filt->scaler_f[1], filt->scaler_f[0] );
		mulVectScalJump_f ( data+oddity, filt->scaler_f[1], len , 2 );		// scalers are inverse!!!
		mulVectScalJump_f ( data+1-oddity, filt->scaler_f[0], len , 2 );
	}

	for (int k=filt->nPredict-1; k>=0; k--)
	{
		XYN_DPRINTF(DEBUG_PROCESS, "%s inverse lifting run %d\n", __FUNCTION__, k );
		// de-update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_f ( data, len, 1-oddity, -1 * filt->coeff_f[1][k] );

		// de-predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_f ( data, len, oddity, -1 * filt->coeff_f[0][k] );
	}
	return OK;
}

/* single lifting stage -
 * if start, first element uses symmetric boundary. */
static int wavLiftingStageInPlace_d ( double* data, long len, int start, register double a )
{
	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, %ld, start? %d) coefficient %f\n",
			__FUNCTION__, data, len, start, a );

	if (start) // do first position d[0] is relevant pass-band, d[0]+= a(d[-1]+d[1])
	{
		*data += a*2*data[1];
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f - start\n", data, *data);
		data+= 2;
		len-=2;
	}
	else		// first position is lowpass, jump over it.
	{
		data++;
		len--;
	}

	while( len > 1 )
	{
		*data += a*(data[-1]+data[1]);
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f (%ld)\n", data, *data, len);
		data += 2;
		len -=2;
	}
	if ( len )		// last position - len was decremented by 2, so if it is ledt 1, do last dot.
	{
		*data += a*2*data[-1];
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f - end\n", data, *data);
	}


	return OK;
}

static int wavLiftingStageInPlace_f ( float* data, long len, int start, register float a )
{
	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, %ld, start? %d) coefficient %f\n",
			__FUNCTION__, data, len, start, a );

	if (start) // do first position d[0] is relevant pass-band, d[0]+= a(d[-1]+d[1])
	{
		*data += a*2*data[1];
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f - start\n", data, *data);
		data+= 2;
		len-=2;
	}
	else		// first position is lowpass, jump over it.
	{
		data++;
		len--;
	}

	while( len > 1 )
	{
		*data += a*(data[-1]+data[1]);
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f (%ld)\n", data, *data, len);
		data += 2;
		len -=2;
	}
	if ( len )		// last position - len was decremented by 2, so if it is ledt 1, do last dot.
	{
		*data += a*2*data[-1];
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f - end\n", data, *data);
	}


	return OK;
}

/* integer lifting transforms */

int liftTransform_i( int* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(len, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) end: (%d %d)\n",
			__FUNCTION__, data, len, oddity, filt->name, sclEnd, wavEnd);

	for (int k=0; k<filt->nPredict; k++)
	{
		// predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_i ( data, len, oddity, filt->coeff_i[0][k] );

		// update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_i ( data, len, 1-oddity, filt->coeff_i[1][k] );
	}

	// packing
	if ( !inPlace )
		wavPack1d_i ( data, len, oddity );

	return OK;
}

int liftInvTransform_i( int* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(len, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) end: (%d %d)\n",
			__FUNCTION__, data, len, oddity, filt->name, sclEnd, wavEnd);

	// unPack
	if ( !inPlace )
		wavUnpack1d_i ( data, len, oddity );

	for (int k=filt->nPredict-1; k>=0; k--)
	{
		XYN_DPRINTF(DEBUG_PROCESS, "%s inverse lifting run %d\n", __FUNCTION__, k );
		// de-update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_i ( data, len, 1-oddity, -1 * filt->coeff_i[1][k] );

		// de-predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_i ( data, len, oddity, -1 * filt->coeff_i[0][k] );
	}
	return OK;
}

static int wavLiftingStageInPlace_i ( int* data, long len, int start, register int a )
{
	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, %ld, start? %d) coefficient %d\n",
			__FUNCTION__, data, len, start, a );

	if (start) // first position d[0] is highpass, d[0]+= a(d[-1]+d[1])
	{
		*data += a*2*data[1] / FixWaveScale;
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4d - start\n", data, *data);
		data+= 2;
		len-=2;
	}
	else		// first position is lowpass, jump over it.
	{
		data++;
		len--;
	}

	while( len > 1 )
	{
		*data += a*(data[-1]+data[1]) / FixWaveScale;
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8d (%ld)\n", data, *data, len);
		data += 2;
		len -=2;
	}
	if ( len )		// last position - len was decremented by 2, so if it is ledt 1, do last dot.
	{
		*data += a*2*data[-1] / FixWaveScale;
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8d - end\n", data, *data);
	}


	return OK;
}


/* rescale */

int rescaleVect_f( float* data, long len )
{
	while (len-- > 0)
		*data++ /= FixWaveScale;
	return OK;
}

int upScaleVect_i( int* data, long len )
{
	while (len-- > 0)
	{
		*data = (*data) << FixWaveShift;
		data ++;
	}
	return OK;
}

int downScaleVect_i( int* data, long len )
{
	while (len-- > 0)
	{
		*data = (*data) >> FixWaveShift;
		data ++;
	}
	return OK;
}


/** vector 1D - floating point version ****************/
int liftTransform_dVect( double* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(vLen, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) end: (%d %d)\n",
			__FUNCTION__, data, vLen, oddity, filt->name, sclEnd, wavEnd);

	/**TODO: check oddity**/
	for (int k=0; k<filt->nPredict; k++)
	{
		// predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_dVect ( data, vLen, vSize, oddity, filt->coeff_f[0][k] );

		// update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"update: \t");
		wavLiftingStageInPlace_dVect ( data, vLen, vSize, 1-oddity, filt->coeff_f[1][k] );
	}

	// packing
	/*if ( !inPlace )
		wavPack1d_fVect ( data, vLen, vSize, oddity );*/

	return OK;
}

int liftTransform_fVect( float* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(vLen, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) end: (%d %d)\n",
			__FUNCTION__, data, vLen, oddity, filt->name, sclEnd, wavEnd);

	/**TODO: check oddity**/
	for (int k=0; k<filt->nPredict; k++)
	{
		// predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_fVect ( data, vLen, vSize, oddity, filt->coeff_f[0][k] );

		// update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"update: \t");
		wavLiftingStageInPlace_fVect ( data, vLen, vSize, 1-oddity, filt->coeff_f[1][k] );
	}

	// packing
	/*if ( !inPlace )
		wavPack1d_fVect ( data, vLen, vSize, oddity );*/

	return OK;
}

int liftInvTransform_fVect( float* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(vLen, oddity);
	int wavEnd = 1-sclEnd;

	// packing
	/*if ( !inPlace )
		wavUnpack1d_fVect ( data, vLen, vSize, oddity );*/

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) end: (%d %d)\n",
			__FUNCTION__, data, vLen, oddity, filt->name, sclEnd,  wavEnd);

	for (int k=filt->nPredict-1; k>=0; k--)
	{
		// update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"de-update: \t");
		wavLiftingStageInPlace_fVect ( data, vLen, vSize, oddity, -1 * filt->coeff_f[1][k] );

		// predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"de-predict: \t");
		wavLiftingStageInPlace_fVect ( data, vLen, vSize, 1-oddity, -1 * filt->coeff_f[0][k] );
	}

	return OK;
}

static int wavLiftingStageInPlace_dVect ( double* data, long vLen, int vSize, int start, register double a )
{
	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, start? %d) coefficient %f\n",
			__FUNCTION__, data, vLen, start, a );

	double *data_p = data, *data_n = data;
	if (start) // first position d[0] is highpass, d[0]+= a(d[-1]+d[1])
	{
		data_n += vSize;
		data_p = data_n;
		addVectorsTriple_d	( data, data_p, data_n, vSize, a, a );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f - start\n", data, *data);
		data   += 2*vSize;
		data_n += 2*vSize;
		vLen -= 2;
	}
	else		// first position is lowpass, jump over it.
	{
		data   += vSize;
		data_n += 2*vSize;
		vLen--;
	}

	while( vLen > 1 )
	{
		addVectorsTriple_d	( data, data_p, data_n, vSize, a, a );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f (%ld)\n", data, *data, vLen);
		data   += 2*vSize;
		data_p += 2*vSize;
		data_n += 2*vSize;
		vLen -= 2;
	}
	if ( vLen )		// last position - len was decremented by 2, so if it is ledt 1, do last dot.
	{
		addVectorsTriple_d	( data, data_p, data_p, vSize, a, a );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f - end\n", data, *data);
	}


	return OK;
}

static int wavLiftingStageInPlace_fVect ( float* data, long vLen, int vSize, int start, register float a )
{
	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, start? %d) coefficient %f\n",
			__FUNCTION__, data, vLen, start, a );

	float *data_p = data, *data_n = data;
	if (start) // first position d[0] is highpass, d[0]+= a(d[-1]+d[1])
	{
		data_n += vSize;
		data_p = data_n;
		addVectorsTriple_f	( data, data_p, data_n, vSize, a, a );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f - start\n", data, *data);
		data   += 2*vSize;
		data_n += 2*vSize;
		vLen -= 2;
	}
	else		// first position is lowpass, jump over it.
	{
		data   += vSize;
		data_n += 2*vSize;
		vLen--;
	}

	while( vLen > 1 )
	{
		addVectorsTriple_f	( data, data_p, data_n, vSize, a, a );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f (%ld)\n", data, *data, vLen);
		data   += 2*vSize;
		data_p += 2*vSize;
		data_n += 2*vSize;
		vLen -= 2;
	}
	if ( vLen )		// last position - len was decremented by 2, so if it is ledt 1, do last dot.
	{
		addVectorsTriple_f	( data, data_p, data_p, vSize, a, a );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8.4f - end\n", data, *data);
	}


	return OK;
}

/** vector 1D - integer lifting ****************/
int liftTransform_iVect( int* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(vLen, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) scl (%d %d) wav (%d %d)\n",
			__FUNCTION__, data, vLen, oddity, filt->name, oddity, sclEnd, 1-oddity, wavEnd);

	for (int k=0; k<filt->nPredict; k++)
	{
		// predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_iVect ( data, vLen, vSize, 1-oddity, filt->coeff_i[0][k] );

		// update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"predict: \t");
		wavLiftingStageInPlace_iVect ( data, vLen, vSize, oddity, filt->coeff_i[1][k] );
	}

	// packing
	/*if ( !inPlace )
		wavPack1d_iVect ( data, vLen, vSize, oddity );*/

	return OK;
}

int liftInvTransform_iVect( int* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	int sclEnd = SCALE_END(vLen, oddity);
	int wavEnd = 1-sclEnd;

	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, oddity %d, %s) end: (%d %d)\n",
			__FUNCTION__, data, vLen, oddity, filt->name, sclEnd, wavEnd);

	// packing
	/*if ( !inPlace )
		wavUnpack1d_iVect ( data, vLen, vSize, oddity );*/

	for (int k=filt->nPredict-1; k>=0; k--)
	{
		// update (low from high)
		XYN_DPRINTF(DEBUG_COEFF,"de-update: \t");
		wavLiftingStageInPlace_iVect ( data, vLen, vSize, oddity, -1 * filt->coeff_i[1][k] );

		// predict (high from low)
		XYN_DPRINTF(DEBUG_COEFF,"de-predict: \t");
		wavLiftingStageInPlace_iVect ( data, vLen, vSize, 1-oddity, -1 * filt->coeff_i[0][k] );
	}

	return OK;
}

static int wavLiftingStageInPlace_iVect ( int* data, long vLen, int vSize, int start, register int a )
{
	XYN_DPRINTF(DEBUG_HEADER, "%s (%p, len %ld, start? %d) coefficient %d\n",
			__FUNCTION__, data, vLen, start, a );

	int *data_p = data, *data_n = data;
	if (start) // first position d[0] is highpass, d[0]+= a(d[-1]+d[1])
	{
		data_n += vSize;
		data_p = data_n;
		addVectorsTriple_i	( data, data_p, data_n, vSize, a, a, FixWaveScale );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8d - start\n", data, *data);
		data   += 2*vSize;
		data_n += 2*vSize;
		vLen -= 2;
	}
	else		// first position is lowpass, jump over it.
	{
		data   += vSize;
		data_n += 2*vSize;
		vLen--;
	}

	while( vLen > 1 )
	{
		addVectorsTriple_i	( data, data_p, data_n, vSize, a, a, FixWaveScale );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8d (%ld)\n", data, *data, vLen);
		data   += 2*vSize;
		data_p += 2*vSize;
		data_n += 2*vSize;
		vLen -= 2;
	}
	if ( vLen )		// last position - len was decremented by 2, so if it is ledt 1, do last dot.
	{
		addVectorsTriple_i	( data, data_p, data_p, vSize, a, a, FixWaveScale );
		XYN_DPRINTF(DEBUG_DETAIL,"%p: %8d - end\n", data, *data);
	}


	return OK;
}

/*** packing/unpacking ***********
 *	used esp. by lifting implementation.
 *	input data, length,
 *	oddity: 1 start with high,
 *			 				--> 1 3 5 7 9... -> start
 *							--> 0 2 4 6 8... -> end
 *			0 start with low,
 *							--> 1 3 5 7 9... -> end
 *							--> 0 2 4 6 8... -> start
 */
static int wavPack1d_d ( double* data, long len, int oddity )
{
	assert ( data!=NULL && len>=0 );
	XYN_DPRINTF(DEBUG_HEADER, "%s Packing\n", __FUNCTION__);

	double* tmp = (double*) malloc(len * sizeof(double));
	memcpy(tmp, data, len * sizeof(double));
	double* tmpCpy = tmp;
	long lowLen = XYN_DIV2(len) - (len % 2) * oddity;		// length of lowpass
	double* outLH[2] = { data, data + lowLen };

	for (long i = 0; i < len; i++, oddity=1-oddity)
		*(outLH[oddity])++ = *tmpCpy++;
	free(tmp);

	return OK;
}

static int wavPack1d_f ( float* data, long len, int oddity )
{
	assert ( data!=NULL && len>=0 );
	XYN_DPRINTF(DEBUG_HEADER, "%s Packing\n", __FUNCTION__);

	float* tmp = (float*) malloc(len * sizeof(float));
	memcpy(tmp, data, len * sizeof(float));
	float* tmpCpy = tmp;
	long lowLen = XYN_DIV2(len) - (len % 2) * oddity;		// length of lowpass
	float* outLH[2] = { data, data + lowLen };

	for (long i = 0; i < len; i++, oddity=1-oddity)
		*(outLH[oddity])++ = *tmpCpy++;
	free(tmp);

	return OK;
}

static int wavUnpack1d_f ( float* data, long len, int oddity )
{
	assert ( data!=NULL && len>=0 );
	XYN_DPRINTF(DEBUG_HEADER, "%s unPacking\n", __FUNCTION__);

	float* tmp = (float*) malloc (len * sizeof(float));
	memcpy (tmp, data, len*sizeof(float));
	long lowLen = XYN_DIV2(len) - (len % 2) * oddity;		// length of lowpass
	float* inpLH[2] = {tmp, tmp + lowLen};

	for (long i=0; i<len; i++, oddity=1-oddity)
		*data++ = *(inpLH[oddity])++;
	free(tmp);

	return OK;
}

static int wavPack1d_i ( int* data, long len, int oddity )
{
	assert ( data!=NULL && len>=0 );
	XYN_DPRINTF(DEBUG_HEADER, "%s Packing\n", __FUNCTION__);

	int* tmp = (int*) malloc(len * sizeof(int));
	memcpy(tmp, data, len * sizeof(int));
	int* tmpCpy = tmp;
	long lowLen = XYN_DIV2(len) - (len % 2) * oddity;		// length of lowpass
	int* outLH[2] = { data, data + lowLen };

	for (long i = 0; i < len; i++, oddity=1-oddity)
		*(outLH[oddity])++ = *tmpCpy++;
	free(tmp);

	return OK;
}

static int wavUnpack1d_i ( int* data, long len, int oddity )
{
	assert ( data!=NULL && len>=0 );
	XYN_DPRINTF(DEBUG_HEADER, "%s unPacking\n", __FUNCTION__);

	int* tmp = (int*) malloc (len * sizeof(int));
	memcpy (tmp, data, len*sizeof(int));
	long lowLen = XYN_DIV2(len) - (len % 2) * oddity;		// length of lowpass
	int* inpLH[2] = {tmp, tmp + lowLen};

	for (long i=0; i<len; i++, oddity=1-oddity)
		*data++ = *(inpLH[oddity])++;
	free(tmp);

	return OK;
}

///////////////////////////Pack 2D//////////////////////////////////////
/*** packing/unpacking ***********
 *	used esp. by lifting implementation.
 *	input data, length,
 *	oddity: 1 start with high,
 *			 				--> 1 3 5 7 9... -> start
 *							--> 0 2 4 6 8... -> end
 *			0 start with low,
 *							--> 1 3 5 7 9... -> end
 *							--> 0 2 4 6 8... -> start
 *
 *	input: data (nC*nR float)
 *	oddity C / R : startHigh on cols / rows
 *	outLLp / outHLp  outLHp outHHp: float** contiguous data nC*nR sized. sizes of specific elements are:
 *
 *			low_nC	high_nC
 *	low_nR	LL  	HL
 *	high_nR	LH  	HH
 *
 *		outLL: lowLenC  * lowLenR
 *		outHL: highLenC * lowLenR
 *		outLH: lowLenC  * highLenR
 *		outHH: highLenC * highLenR
 *
 *		note that only *outLLp need to be freed, since the output buffers are contiguous.
 */
int wavPack_lineBuffers_f ( float* data, long nC, long nR, int oddityC, bool oddityR, float** outLLp, float** outHLp, float** outLHp, float** outHHp, long *low_nC, long *low_nR, long *high_nC, long *high_nR )
{
	assert ( data!=NULL && nC>=0 && nR>=0 );
	XYN_DPRINTF(DEBUG_HEADER, "%s Packing\n", __FUNCTION__);

	long lowLenR = XYN_DIV2(nR) - (nR % 2) * oddityR;		// lowpass nRows
	long highLenR = nR-lowLenR;

	long lowLenC = XYN_DIV2(nC)  - (nC % 2) * oddityC;		// lowpass nCols
	long highLenC = nC-lowLenC;

	//float* outLH[2] = { data, data + lowLen };
	//float* outLL = (float*) calloc(lowLenC  * lowLenR  , sizeof(float));
	//float* outLH = (float*) calloc(highLenC * lowLenR  , sizeof(float));
	//float* outHL = (float*) calloc(lowLenC  * highLenR , sizeof(float));
	//float* outHH = (float*) calloc(highLenC * highLenR , sizeof(float));

	// allocate contiguous buffer, and move band pass elements relative to it.
	float* outLL = (float*) calloc(nC * nR, sizeof(float));		// save lowLenC  * lowLenR floats
	float* outHL = outLL + lowLenC  * lowLenR;					// save highLenC * lowLenR floats
	float* outLH = outHL + highLenC * lowLenR;					// save lowLenC  * highLenR floats
	float* outHH = outLH + lowLenC  * highLenR;					// has  highLenC * highLenR floats

	float* outLineL[2] = { outLL, outLH };
	float* outLineH[2] = { outHL, outHH };

	XYN_DPRINTF(DEBUG_HEADER, "\n sizes: low nC nR high nC nR \n %ld %ld %ld %ld \n pointers: %p %p %p %p\n", lowLenC, lowLenR, highLenC, highLenR, outLL, outHL, outLH, outHH);
	//displayArrayNM_f(data, nC, nR, stdout );
	XYN_DPRINTF(DEBUG_HEADER, "\n------------------LL-----------------------\n");
	//memset(outLL, 1, lowLenC  * lowLenR  * sizeof(float));
	//displayArrayNM_f(outLL, lowLenC, lowLenR, stdout );

	for (long i = 0; i < nR; i++, oddityR = 1-oddityR)
	{
		memcpy(outLineL[oddityR], data, 		  lowLenC *sizeof(float));
		memcpy(outLineH[oddityR], data + lowLenC, highLenC*sizeof(float));


		XYN_DPRINTF(DEBUG_HEADER, "------------------Line %ld: (%d) %p %p %p %p (%ld %ld)\n", i, oddityR, outLineL[0], outLineH[0], outLL, outHL, outLineL[0]-outLL, outLineH[0]-outHL);
		XYN_DPRINTF(DEBUG_HEADER, "------------------Line %ld: (%d) %p %p %p %p (%ld %ld)\n", i, oddityR, outLineL[1], outLineH[1], outLH, outHH, outLineL[1]-outLH, outLineH[1]-outHH);
		//displayArray_f(data, nC,  stdout );
		//displayArray_f(data + lowLenC, highLenC,  stdout );
		//displayArray_f(outLineL[oddityR], lowLenC,  stdout );
		//displayArray_f(outLineH[oddityR], highLenC, stdout );
		XYN_DPRINTF(DEBUG_HEADER, "\n");

		outLineL[oddityR] += lowLenC ;
		outLineH[oddityR] += highLenC;
		data += nC;
		//*(outLH[oddity])++ = *tmpCpy++;
	}

	*outLLp = outLL;
	*outLHp = outLH;
	*outHLp = outHL;
	*outHHp = outHH;

	*low_nC = lowLenC;
	*low_nR = lowLenR;
	*high_nC = highLenC;
	*high_nR = highLenR;
	XYN_DPRINTF(DEBUG_HEADER, "\n sizes: low nC nR high nC nR \n %ld %ld %ld %ld \n pointers: %p %p %p %p\n", *low_nC, *low_nR, *high_nC, *high_nR, outLL, outHL, outLH, outHH);

/*
	printf("\n\n -------- LL ---------------\n");
	displayArrayNM_f(outLL, *low_nC, *low_nR, stdout );

	printf("\n\n -------- LH ---------------\n");
	displayArrayNM_f(outLH, *high_nC, *low_nR, stdout );

	printf("\n\n -------- HL ---------------\n");
	displayArrayNM_f(outHL, *low_nC, *high_nR, stdout );

	printf("\n\n -------- HH ---------------\n");
	displayArrayNM_f(outHH, *high_nC, *high_nR, stdout );
*/

	return OK;
}

int wavPack_lineBuffers_d ( double* data, long nC, long nR, int oddityC, bool oddityR, double** outLLp, double** outHLp, double** outLHp, double** outHHp, long *low_nC, long *low_nR, long *high_nC, long *high_nR )
{
	assert ( data!=NULL && nC>=0 && nR>=0 );
	XYN_DPRINTF(DEBUG_HEADER, "%s Packing\n", __FUNCTION__);

	long lowLenR = XYN_DIV2(nR) - (nR % 2) * oddityR;		// lowpass nRows
	long highLenR = nR-lowLenR;

	long lowLenC = XYN_DIV2(nC)  - (nC % 2) * oddityC;		// lowpass nCols
	long highLenC = nC-lowLenC;

	//float* outLH[2] = { data, data + lowLen };
	//float* outLL = (float*) calloc(lowLenC  * lowLenR  , sizeof(float));
	//float* outLH = (float*) calloc(highLenC * lowLenR  , sizeof(float));
	//float* outHL = (float*) calloc(lowLenC  * highLenR , sizeof(float));
	//float* outHH = (float*) calloc(highLenC * highLenR , sizeof(float));

	// allocate contiguous buffer, and move band pass elements relative to it.
	double* outLL = (double*) calloc(nC * nR, sizeof(double));		// save lowLenC  * lowLenR floats
	double* outHL = outLL + lowLenC  * lowLenR;					// save highLenC * lowLenR floats
	double* outLH = outHL + highLenC * lowLenR;					// save lowLenC  * highLenR floats
	double* outHH = outLH + lowLenC  * highLenR;					// has  highLenC * highLenR floats

	double* outLineL[2] = { outLL, outLH };
	double* outLineH[2] = { outHL, outHH };

	XYN_DPRINTF(DEBUG_HEADER, "\n sizes: low nC nR high nC nR \n %ld %ld %ld %ld \n pointers: %p %p %p %p\n", lowLenC, lowLenR, highLenC, highLenR, outLL, outHL, outLH, outHH);
	//displayArrayNM_f(data, nC, nR, stdout );
	XYN_DPRINTF(DEBUG_HEADER, "\n------------------LL-----------------------\n");
	//memset(outLL, 1, lowLenC  * lowLenR  * sizeof(float));
	//displayArrayNM_f(outLL, lowLenC, lowLenR, stdout );

	for (long i = 0; i < nR; i++, oddityR = 1-oddityR)
	{
		memcpy(outLineL[oddityR], data, 		  lowLenC *sizeof(double));
		memcpy(outLineH[oddityR], data + lowLenC, highLenC*sizeof(double));


		XYN_DPRINTF(DEBUG_HEADER, "------------------Line %ld: (%d) %p %p %p %p (%ld %ld)\n", i, oddityR, outLineL[0], outLineH[0], outLL, outHL, outLineL[0]-outLL, outLineH[0]-outHL);
		XYN_DPRINTF(DEBUG_HEADER, "------------------Line %ld: (%d) %p %p %p %p (%ld %ld)\n", i, oddityR, outLineL[1], outLineH[1], outLH, outHH, outLineL[1]-outLH, outLineH[1]-outHH);
		//displayArray_f(data, nC,  stdout );
		//displayArray_f(data + lowLenC, highLenC,  stdout );
		//displayArray_f(outLineL[oddityR], lowLenC,  stdout );
		//displayArray_f(outLineH[oddityR], highLenC, stdout );
		XYN_DPRINTF(DEBUG_HEADER, "\n");

		outLineL[oddityR] += lowLenC ;
		outLineH[oddityR] += highLenC;
		data += nC;
		//*(outLH[oddity])++ = *tmpCpy++;
	}

	*outLLp = outLL;
	*outLHp = outLH;
	*outHLp = outHL;
	*outHHp = outHH;

	*low_nC = lowLenC;
	*low_nR = lowLenR;
	*high_nC = highLenC;
	*high_nR = highLenR;
	XYN_DPRINTF(DEBUG_HEADER, "\n sizes: low nC nR high nC nR \n %ld %ld %ld %ld \n pointers: %p %p %p %p\n", *low_nC, *low_nR, *high_nC, *high_nR, outLL, outHL, outLH, outHH);

/*
	printf("\n\n -------- LL ---------------\n");
	displayArrayNM_f(outLL, *low_nC, *low_nR, stdout );

	printf("\n\n -------- LH ---------------\n");
	displayArrayNM_f(outLH, *high_nC, *low_nR, stdout );

	printf("\n\n -------- HL ---------------\n");
	displayArrayNM_f(outHL, *low_nC, *high_nR, stdout );

	printf("\n\n -------- HH ---------------\n");
	displayArrayNM_f(outHH, *high_nC, *high_nR, stdout );
*/

	return OK;
}

/** 2D  ********************************************************************************/

int liftTransform2d_d ( double* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	assert( nC>0 && nR>0 && data!=NULL && filt!=NULL );
	double* dataCpy = data;
	for ( int j=0; j<nR; j++ )
	{
		liftTransform_d ( dataCpy, nC, oddityC, inPlace, noRescale, filt );
		dataCpy += nC;
	}
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done columns 2d\n", __FUNCTION__);
	//displayArrayNM_f ( data, nC, nR, stdout );

	liftTransform_dVect( data, nR, nC, oddityR, inPlace, noRescale, filt );
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done rows 2d\n", __FUNCTION__);
	//displayArrayNM_f ( data, nC, nR, stdout );

	return OK;
}

int liftTransform2d_f ( float* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	assert( nC>0 && nR>0 && data!=NULL && filt!=NULL );
	float* dataCpy = data;
	for ( int j=0; j<nR; j++ )
	{
		liftTransform_f ( dataCpy, nC, oddityC, inPlace, noRescale, filt );
		dataCpy += nC;
	}
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done columns 2d\n", __FUNCTION__);
	//displayArrayNM_f ( data, nC, nR, stdout );

	liftTransform_fVect( data, nR, nC, oddityR, inPlace, noRescale, filt );
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done rows 2d\n", __FUNCTION__);
	//displayArrayNM_f ( data, nC, nR, stdout );

	return OK;
}

int liftInvTransform2d_f ( float* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	assert( nC>0 && nR>0 && data!=NULL && filt!=NULL );
	float* dataCpy = data;
	for ( int j=0; j<nR; j++ )
	{
		liftInvTransform_f ( dataCpy, nC, oddityC, inPlace, noRescale, filt );
		dataCpy += nC;
	}
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done columns 2d\n", __FUNCTION__);
	//displayArrayNM_f ( data, nC, nR, stdout );

	liftInvTransform_fVect( data, nR, nC, oddityR, inPlace, noRescale, filt );
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done rows 2d\n", __FUNCTION__);
	//displayArrayNM_f ( data, nC, nR, stdout );

	return OK;
}

int liftTransform2d_i ( int* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	assert( nC>0 && nR>0 && data!=NULL && filt!=NULL );
	int* dataCpy = data;
	for ( int j=0; j<nR; j++ )
	{
		liftTransform_i ( dataCpy, nC, oddityC, inPlace, noRescale, filt );
		dataCpy += nC;
	}
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done columns 2d\n", __FUNCTION__);
	//displayArrayNM_i ( data, nC, nR, stdout );

	liftTransform_iVect( data, nR, nC, oddityR, inPlace, noRescale, filt );
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done rows 2d\n", __FUNCTION__);
	//displayArrayNM_i ( data, nC, nR, stdout );

	return OK;
}

int liftInvTransform2d_i ( int* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt )
{
	assert( nC>0 && nR>0 && data!=NULL && filt!=NULL );
	int* dataCpy = data;
	for ( int j=0; j<nR; j++ )
	{
		liftInvTransform_i ( dataCpy, nC, oddityC, inPlace, noRescale, filt );
		dataCpy += nC;
	}
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done columns 2d\n", __FUNCTION__);
	//displayArrayNM_i ( data, nC, nR, stdout );

	liftInvTransform_iVect( data, nR, nC, oddityR, inPlace, noRescale, filt );
	XYN_DPRINTF(DEBUG_PROCESS, "%s: done rows 2d\n", __FUNCTION__);
	//displayArrayNM_i ( data, nC, nR, stdout );

	return OK;
}


/***** static fuctions **********/
/* initialize filters */
static void wavInitHaar()
{
	waveFilterPtr wfp = &wavHaar;

	wfp->name = "Haar";
	wfp->len_anal[0] = wfp->len_anal[1] = wfp->len_synt[0] = wfp->len_synt[1] = 2;
	wfp->pos_anal[0] = wfp->pos_anal[1] = wfp->pos_synt[0] = wfp->pos_synt[1] = 0;

	wfp->f_anal[0][0] = wfp->f_anal[0][1] = wfp->f_anal[1][0] = wfp->f_anal[1][1] = M_SQRT1_2;
	wfp->f_synt[0][0] = wfp->f_synt[0][1] = wfp->f_synt[1][0] = wfp->f_synt[1][1] = M_SQRT1_2;

	wavCompleteWavFilterData ( wfp );

	return;
}

static void wavInitTest()
{
	waveFilterPtr wfp = &wavTest;
	wfp->name = "Test";

	wfp->len_anal[0] = wfp->len_synt[1] = 5;
	wfp->len_synt[0] = wfp->len_anal[1] = 7;
	wfp->pos_anal[0] = wfp->pos_synt[1] = 2;
	wfp->pos_synt[0] = wfp->pos_anal[1] = 3;

	wfp->f_anal[0][2] =  1;
	wfp->f_anal[0][1] =  2.2;
	wfp->f_anal[0][0] =  3.03;

	wfp->f_anal[1][3] =  1;
	wfp->f_anal[1][2] =  3.03;
	wfp->f_anal[1][1] =  4.004;
	wfp->f_anal[1][0] =  5.0005;

	wavCompleteWavFilterData ( wfp );

	return;
}

static void wavInitCDF53()
{
	waveFilterPtr wfp = &wavCDF53;
	wfp->name = "CDF-5-3";

	wfp->len_anal[0] = wfp->len_synt[1] = 5;
	wfp->len_synt[0] = wfp->len_anal[1] = 3;
	wfp->pos_anal[0] = wfp->pos_synt[1] = 2;
	wfp->pos_synt[0] = wfp->pos_anal[1] = 1;

	wfp->f_anal[0][0] = -0.25 * M_SQRT1_2;
	wfp->f_anal[0][1] =  0.5  * M_SQRT1_2;
	wfp->f_anal[0][2] =  1.5  * M_SQRT1_2;

	wfp->f_anal[1][0] = -0.5 * M_SQRT1_2;
	wfp->f_anal[1][1] =  1.0 * M_SQRT1_2;

	wavCompleteWavFilterData ( wfp );

	return;
}

static void wavInitCDF97()
{
	waveFilterPtr wfp = &wavCDF97;
	wfp->name = "CDF-9-7";

	wfp->len_anal[0] = wfp->len_synt[1] = 9;
	wfp->len_synt[0] = wfp->len_anal[1] = 7;
	wfp->pos_anal[0] = wfp->pos_synt[1] = 4;
	wfp->pos_synt[0] = wfp->pos_anal[1] = 3;
 /*filters:
  * analysis low
  *  0.85269865321930, 0.37740268810913, -0.11062402748951, -0.02384929751586,  0.03782879857992
  * high
  *  0.78848487220618, 0.41809244072573, -0.04068975261660, -0.06453905013246
  */
	wfp->f_anal[0][0] =  0.03782879857992;
	wfp->f_anal[0][1] = -0.02384929751586;
	wfp->f_anal[0][2] = -0.11062402748951;
	wfp->f_anal[0][3] =  0.37740268810913;
	wfp->f_anal[0][4] =  0.85269865321930;

	wfp->f_anal[1][0] =  0.06453905013246;
	wfp->f_anal[1][1] = -0.04068975261660;
	wfp->f_anal[1][2] = -0.41809244072573;
	wfp->f_anal[1][3] =  0.78848487220618;

	wavCompleteWavFilterData ( wfp );

	//wavGenFilterFold(ftp);
	return;
}



static void wavInitLiftCDF53()
{
	waveLiftPtr wlp = &wavLiftCDF53;
	wlp->name = "CDF-5-3";

	wlp->nPredict = 1;				//1 pass
	wlp->coeff_f[0][0] = -0.5;		//predict
	wlp->coeff_f[1][0] = 0.25;		//update

	wlp->scaler_f[0] = M_SQRT2;		//lowpass scale
	wlp->scaler_f[1] = M_SQRT1_2;	//highpass scale

	/* prepare integer coefficients */
	wavCoeffsIntFromFloat ( wlp->coeff_i[0], wlp->coeff_f[0], wlp->nPredict );
	wavCoeffsIntFromFloat ( wlp->coeff_i[1], wlp->coeff_f[1], wlp->nPredict );
	return;
}


static void wavInitLiftCDF97()
{
	/* coefficients are -1.586134342, -0.05298011854, 0.8829110762, 0.4435068522 */
	waveLiftPtr wlp = &wavLiftCDF97;
	wlp->name = "CDF-9-7";

	wlp->nPredict = 2;
	wlp->coeff_f[0][0] = -1.586134342;
	wlp->coeff_f[1][0] = -0.05298011854;
	wlp->coeff_f[0][1] =  0.8829110762;
	wlp->coeff_f[1][1] =  0.4435068522;

	wlp->scaler_f[0] = 1.149604398;
	wlp->scaler_f[1] = 1.0 / wlp->scaler_f[0];

	wavCoeffsIntFromFloat ( wlp->coeff_i[0], wlp->coeff_f[0], wlp->nPredict );
	wavCoeffsIntFromFloat ( wlp->coeff_i[1], wlp->coeff_f[1], wlp->nPredict );

	return;
}


static void genPolyFilterFromWave ( waveFilterPtr wfp, wavePolyPtr pwp )
{
	/*
	for (int i=0; i<=1; i++)
	{
		pwp->pos_anal[i][0] = DIV2( wfp->pos_anal[i] );
		pwp->pos_anal[i][1] = ( wfp->pos_anal[i] ) / 2;

		pwp->pos_synt[i][0] = DIV2( wfp->pos_synt[i] );
		pwp->pos_synt[i][1] = ( wfp->pos_synt[i] ) / 2;

		pwp->len_anal[i][0] = DIV2( wfp->len_anal[i] );
		pwp->len_anal[i][1] = ( wfp->len_anal[i] ) / 2;

		pwp->len_synt[i][0] = DIV2( wfp->len_synt[i] );
		pwp->len_synt[i][1] = ( wfp->len_synt[i] ) / 2;
	}
	*/
	for (int i=0; i<=1; i++)				// i=0 lowpass; i=1 highpass
	for (int k=0; k<=1; k++)				// k=0 same oddity as zero (even places); k=1 odd places
	{
		bool m_anal = i ^ (wfp->pos_anal[k] %2);		// relative movement of PP filter is zPos xor i
		bool m_synt = i ^ (wfp->pos_synt[k] %2);

		pwp->pos_anal[i][k] = ( wfp->pos_anal[i] + 1-m_anal) / 2;
		pwp->pos_synt[i][k] = ( wfp->pos_synt[i] + 1-m_synt) / 2;
		pwp->len_anal[i][k] = ( wfp->len_anal[i] + 1-m_anal) / 2;
		pwp->len_synt[i][k] = ( wfp->len_synt[i] + 1-m_synt) / 2;

		XYN_DPRINTF(DEBUG_PROCESS, "\n%s for positions %d: positions change by %d/%d\n"
				"positions - anal/synt %d(%d) / %d(%d)\n"
				"lengths   - anal/synt %d(%d) / %d(%d)\n",
				__FUNCTION__, i, m_anal, m_synt,
				pwp->pos_anal[i][k], wfp->pos_anal[k], pwp->pos_synt[i][k], wfp->pos_synt[k],
				pwp->len_anal[i][k], wfp->len_anal[k], pwp->len_synt[i][k], wfp->len_synt[k] );

		cpyJump_f( wfp->f_anal[i] + m_anal, 2, pwp->f_anal[i][k], 1, pwp->len_anal[i][k] );
		XYN_DPRINTF(DEBUG_DETAIL, "anal[%d][%d]\t", i, k); displayArray_f(pwp->f_anal[i][k], pwp->len_anal[i][k], stdout);
		cpyJump_f( wfp->f_synt[i] + m_synt, 2, pwp->f_synt[i][k], 1, pwp->len_synt[i][k] );
		XYN_DPRINTF(DEBUG_DETAIL, "synt[%d][%d]\t", i, k); displayArray_f(pwp->f_synt[i][k], pwp->len_synt[i][k], stdout);

		for( int j=0; j<wfp->len_anal[i]; j++ )
		{
			pwp->len_anal_folded[i][k][j] = ( wfp->len_anal_folded[i][j] ) / 2;

			XYN_DPRINTF(DEBUG_PROCESS, "\tfold positions - anal %d %d\n\tfold lengths   - anal %d %d\n",
							pwp->pos_anal[i][0], pwp->pos_anal[i][1], pwp->len_anal_folded[i][0][j], pwp->len_anal_folded[i][1][j] );

			displayArray_f(wfp->f_anal_folded[i][j], wfp->len_anal_folded[i][j], stdout);
			cpyJump_f( wfp->f_anal_folded[i][j] + m_anal, 2, pwp->f_anal_folded[i][k][j], 1, pwp->len_anal_folded[i][k][j] );
		}
		for( int j=0; j<wfp->len_synt[i]; j++ )
			{
				pwp->len_synt_folded[i][k][j] = ( wfp->len_synt_folded[i][j] ) / 2;

				XYN_DPRINTF(DEBUG_PROCESS, "\tfold positions - synt %d %d\n\tfold lengths   - synt %d %d\n",
								pwp->pos_synt[i][0], pwp->pos_synt[i][1], pwp->len_synt_folded[i][0][j], pwp->len_synt_folded[i][1][j]);

				displayArray_f(wfp->f_synt_folded[i][j], wfp->len_synt_folded[i][j], stdout);
				cpyJump_f( wfp->f_synt_folded[i][j] + m_synt, 2, pwp->f_synt_folded[i][k][j], 1, pwp->len_synt_folded[i][k][j] );
			}
	}


	return;
}


/** filter generators: symmetric copy, QMF, etc... **/

/** QMF generator: change sign on every other coefficient **/
static int wavFilterQMF_i(int* out, int* inp, int len, int zeroPos)
{
	assert (inp && out );
	int QMFSign[2] = {1,-1};
	int oddity = zeroPos%2;

	while (len-- > 0)
	{
		*out++ = QMFSign[oddity] * (*inp++);
		oddity = 1-oddity;
	}
	return OK;
}

static int wavFilterQMF_f(float* out, float* inp, int len, int zeroPos)
{
	assert (inp && out );
	int QMFSign[2] = {1,-1};
	int oddity = zeroPos%2;

	while (len-- > 0)
	{
		*out++ = QMFSign[oddity] * (*inp++);
		oddity = 1-oddity;
	}
	return OK;
}

/*
 * build symmetric folds - 0 to len-1:
 * suppose filter is  -5 -4 -3 -2 -1 0 1 2 3 4 5
 * then folds are at +- 1 2 3 4 5 (fold 0 is the initial filter)
 * fold at 1 is: -5 -4 -3 -2 -1 0 1 2 3+5 4
 * fold at 2 is: -5 -4 -3 -2 -1 0 1+5 2+4 3
 * 12345 zero at -2, position 0: if (zero+pos)>=0 do nothing.
 *
 */
static int wavFilterFold_f(float (*out)[FILTER_LEN], int* outLen, float* inp, int len, int zeroPos)
{
	assert (inp && out && FILTER_LEN>=len );

	int i=0, k_i = zeroPos;
	XYN_DPRINTF(DEBUG_DETAIL,"left-side folding - %d to %d of %d\n", zeroPos, 0, len);
	for (; k_i>=0; i++, k_i--)
	{
		memcpy(out[i], inp+k_i, (len-k_i)*sizeof(float));
		//displayArray_f(out[i], len, stdout);
		for (int j=1, k_j=1; k_j<=k_i; k_j++, j++ )
			out[i][j] += inp[k_i-k_j];

		outLen[i] = len - k_i;//zeroPos + i + 1;
		//displayArray_f(out[i], len, stdout);NL;
	}

	XYN_DPRINTF(DEBUG_DETAIL,"right-side folding - i %d k_i %d of %d\n", i, k_i, len);
	for (; i<len; i++, k_i--)
	{
		outLen[i] = len+k_i;

		memcpy(out[i], inp, (len+k_i)*sizeof(float));
		//displayArray_f(out[i], len, stdout);
		for (int j=outLen[i]-2, k_j=len+k_i; k_j<len; j--, k_j++ )
			out[i][j] += inp[k_j];
		//displayArray_f(out[i], len, stdout);NL;


	}

	return OK;
}

static int wavFilterSym_f(float* data, int len, int zeroPos)
{
	assert (data );
	XYN_DPRINTF(DEBUG_HEADER,"%s generating symmetric filter complements from (-1)->(%d) (len %d)\n", __FUNCTION__, zeroPos, len);

	for (int i=0; i<=zeroPos; i++)
		data[len-1-i] = data[i];

	return OK;
}

static int wavCoeffsIntFromFloat ( int* int_coeffs, float* float_coeffs, int n )
{
	for ( int i=0; i<n; i++ )
		int_coeffs[i] = (int) ( float_coeffs[i] * FixWaveScale );
	return OK;
}

static void wavCompleteWavFilterData ( waveFilterPtr wfp )
{
	for (int i=0; i<=1; i++)	// even/odd (low/high)
	{
		wavFilterSym_f(wfp->f_anal[i], wfp->len_anal[i], wfp->pos_anal[i] );
		wavFilterQMF_f(wfp->f_synt[1-i], wfp->f_anal[i], wfp->len_anal[i], wfp->pos_anal[i] );
	}
	for (int i=0; i<=1; i++)	// even/odd (low/high)
	{
		wavCoeffsIntFromFloat ( wfp->i_anal[i], wfp->f_anal[i], wfp->len_anal[i] );
		wavCoeffsIntFromFloat ( wfp->i_synt[i], wfp->f_synt[i], wfp->len_synt[i] );

		XYN_DPRINTF(DEBUG_DETAIL, "folding f_anal[%d]. len %d. zeroPos %d\n", i, wfp->len_anal[i], wfp->pos_anal[i]);
//			displayArray_f(wfp->f_anal[i], wfp->len_anal[i], stdout);

		wavFilterFold_f	( wfp->f_anal_folded[i], wfp->len_anal_folded[i], wfp->f_anal[i], wfp->len_anal[i], wfp->pos_anal[i] );

		XYN_DPRINTF(DEBUG_DETAIL, "folding f_synt[%d]. len %d. zeroPos %d\n", i, wfp->len_synt[i], wfp->pos_synt[i]);
//			displayArray_f(wfp->f_synt[i], wfp->len_synt[i], stdout);
		wavFilterFold_f	( wfp->f_synt_folded[i], wfp->len_synt_folded[i], wfp->f_synt[i], wfp->len_synt[i], wfp->pos_synt[i] );

		for (int j=0; j<wfp->len_anal[i]; j++)
			wavCoeffsIntFromFloat ( wfp->i_anal_folded[i][j], wfp->f_anal_folded[i][j], wfp->len_anal[i] );
		for (int j=0; j<wfp->len_synt[i]; j++)
			wavCoeffsIntFromFloat ( wfp->i_synt_folded[i][j], wfp->f_synt_folded[i][j], wfp->len_synt[i] );
	}

	return;
}













/* end */
