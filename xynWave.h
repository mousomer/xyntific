/*
 * xynWave.h
 *
 *  Created on: Feb 5, 2011
 *      Author: mousomer
 */

#ifndef XYNWAVE_H_
#define XYNWAVE_H_

enum waveType { xynWaveleType=0, xynLiftingType=1, xynPolyPhaseType=2 };

struct waveletFilter{
	char* name;
	/* integer filter */
	int i_anal[2][FILTER_LEN];			/* lowpass/highpass analysis  filter */
	int i_synt[2][FILTER_LEN];			/* lowpass/highpass synthesis filter */
	int i_anal_folded[2][FILTER_LEN][FILTER_LEN];
	int i_synt_folded[2][FILTER_LEN][FILTER_LEN];
	/* floating point filter */
	float f_anal[2][FILTER_LEN];		/* lowpass/highpass analysis  filter */
	float f_synt[2][FILTER_LEN];		/* lowpass/highpass synthesis filter */
	float f_anal_folded[2][FILTER_LEN][FILTER_LEN];
	float f_synt_folded[2][FILTER_LEN][FILTER_LEN];

	/* filter properties - no need for boundary zeroPosition (it is 0 for lower boundaries, normal for upper boundaries) */
	int len_anal[2], len_synt[2];	// filter lengths
	int pos_anal[2], pos_synt[2];	// starting filter positions relative to 0
	int len_anal_folded[2][FILTER_LEN], len_synt_folded[2][FILTER_LEN];
};

struct wavePolyFilter{
	char* name;
	/* floating point filters */
	float f_anal[2][2][POLY_LEN];	// lowpass/highpass even/odd - analysis polyphase floating pofloat filters
	float f_synt[2][2][POLY_LEN];	// lowpass/highpass even/odd - synthesis polyphase floating pofloat filters
	float f_anal_folded[2][2][POLY_LEN][FILTER_LEN];	// folded boundary filters - analysis
	float f_synt_folded[2][2][POLY_LEN][FILTER_LEN];	// folded boundary filters - synthesis

	/* integer filters */
	int i_anal[2][2][POLY_LEN];	// lowpass/highpass even/odd  analysis polyphase integer filters
	int i_synt[2][2][POLY_LEN];	// lowpass/highpass even/odd  analysis polyphase integer filters
	int i_anal_folded[2][2][POLY_LEN][FILTER_LEN];	// folded boundary filters - analysis
	int i_synt_folded[2][2][POLY_LEN][FILTER_LEN];	// folded boundary filters - synthesis

	/* filter properties */
	int len_anal[2][2], len_synt[2][2];	// filter lengths
	int pos_anal[2][2], pos_synt[2][2];	// starting filter positions relative to 0
	int len_anal_folded[2][2][FILTER_LEN], len_synt_folded[2][2][FILTER_LEN];
};

struct waveLiftFilter{
	char* name;
	int nPredict;	// prediction level

	int 	scaler_i[2];			// integer scaling coefficients
	int		coeff_i[2][POLY_LEN];	// integer prediction/update coefficients
	float 	scaler_f[2]; 			// float   scaling coefficients
	float 	coeff_f[2][POLY_LEN];	// float   prediction/update coefficients
};

struct dwtInfo{
	waveFilterPtr fwp;
	wavePolyPtr	  pwp;
	waveLiftPtr	  lwp;
	enum waveType 	wType;

	bool oddity, inPlace;
	int jump;
};


waveLiftPtr		wavGetLiftWaveName 	( const char* name );
wavePolyPtr		wavGetPolyWaveName	( const char* name );
waveFilterPtr	wavGetWaveletName	( const char* name );

dwtInfPtr 		wavGenDWTInfo 		( const char* wavName, enum waveType wtype, bool oddity, bool inPlace );
int wavPrint ( waveFilterPtr filt, FILE* fp );
int wavLiftPrint ( waveLiftPtr filt, FILE* fp );
int wavPolyPrint ( wavePolyPtr filt, FILE* fp );


/** integer transforms **/
int wavTransform_i			( int* dataInp, int* dataOut, long len, bool oddity, waveFilterPtr filt );
int wavInvTransform_i		( int* dataInp, int* dataOut, long len, bool oddity, waveFilterPtr filt );
int wavTransformInPlace_i	( int* data, long len, bool oddity, waveFilterPtr filt );
int wavInvTransformInPlace_i( int* data, long len, bool oddity, waveFilterPtr filt );
int liftTransform_i			( int* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftInvTransform_i		( int* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );

int downScaleVect_i			( int* data, long len );
int upScaleVect_i			( int* data, long len );

int liftTransform_iVect		( int* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftInvTransform_iVect	( int* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftTransform2d_i 		( int* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftInvTransform2d_i	( int* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt );

/** float transforms **/
int wavTransform_f			( float* dataInp, float* dataOut, long len, bool oddity, bool inPlace, waveFilterPtr filt );
int wavInvTransform_f		( float* dataInp, float* dataOut, long len, bool oddity, bool inPlace, waveFilterPtr filt );

//int wavTransformInPlace_f	( float* dataInp, float* dataOut, long len, bool oddity, waveFilterPtr filt );
//int wavInvTransformInPlace_f( float* dataInp, float* dataOut, long len, bool oddity, waveFilterPtr filt );
/*int wavInvTransformHiLow_f	( float* dataOut, long len, float* inpLow, long lenLow,
							  float* inpHigh, long lenHigh, bool oddity, waveFilterPtr filt );
int wavTransformHiLow_f		( float* dataInp, long len, float* outLow, long lenLow,
							float* outHigh, long lenHigh, bool oddity, waveFilterPtr filt );
*/
int liftTransform_d			( double* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftTransform_f			( float* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftInvTransform_f		( float* data, long len, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );

int rescaleVect_f			( float* data, long len );

int liftTransform_fVect		( float* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftTransform_dVect		( double* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftInvTransform_fVect	( float* data, long vLen, int vSize, bool oddity, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftTransform2d_f 		( float* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftTransform2d_d 		( double* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt );
int liftInvTransform2d_f	( float* data, int nC, int nR, bool oddityC, bool oddityR, bool inPlace, bool noRescale, waveLiftPtr filt );
//int liftTransform_f		( float* dataInp, float* dataOut, long len, bool oddity, waveLiftPtr filt );
//int liftInvTransform_f	( float* dataInp, float* dataOut, long len, bool oddity, waveLiftPtr filt );



/** packing / unpacking **/
int wavPack1d_fVect		( float* data, long vLen, int vSize, int oddity );
int wavUnpack1d_fVect 	( float* data, long vLen, int vSize, int oddity );
int wavPack1d_iVect		( int* data, long vLen, int vSize, int oddity );
int wavUnpack1d_iVect 	( int* data, long vLen, int vSize, int oddity );

/** 2d line packing - mallocing the line buffers **/
int wavPack_lineBuffers_f ( float* data, long nC, long nR, int oddityC, bool oddityR, float** outLLp, float** outHLp, float** outLHp, float** outHHp, long *low_nC, long *low_nR, long *high_nC, long *high_nR );
int wavPack_lineBuffers_d ( double* data, long nC, long nR, int oddityC, bool oddityR, double** outLLp, double** outHLp, double** outLHp, double** outHHp, long *low_nC, long *low_nR, long *high_nC, long *high_nR );






#endif /* XYNWAVE_H_ */
