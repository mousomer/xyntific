/*
 * dwt2d.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#include "libxyn.h"

#define USG_STRING "[-v %:] [-out %s:outName] [-noScl %:] [-startHigh %:] [-nC %d:nC] %d:nR %s:waveName %s:inpName"
xynString waveName, inpName, outName, outNameLL, outNameLH, outNameHL, outNameHH;
int verbose=0, startHigh=0, noRescale=0, nC=0, nR=0;

#define PREPARE_DATA(_t_) {		\
		XYN_DPRINTF(DEBUG_PROCESS,"\npreparing data. Length %d/=%d\n", len, nBytes);	\
		_t_##Input  = (_t_*) malloc(nBytes * len);	\
		_t_##Output = (_t_*) malloc(nBytes * len);	\
		assert(_t_##Input && _t_##Output );			\
		fread(_t_##Input, nBytes, len, fpInp);		\
		memcpy(_t_##Output, _t_##Input, len*nBytes );	}


int main(int argc, char* argv[])
{
	if (parseParameters(argc, argv, USG_STRING,
			&verbose, outName, &noRescale, &startHigh, &nC, &nR, waveName, inpName))
		return (FAIL);

	/* transform types:
			0: shuffled,	non lift
			1: inPlace, 	non lift
			2: shuffled, 	lifting implementation
			3: inPlace, 	lifting implementation	*/
	//int transformType = 4*(dataType[0]=='f') + 2*(lifting>0) + (inPlace>0);

	long len = getFileSize(inpName);
	int nBytes = sizeof(float);
	len /= nBytes;
	if (nC == 0)
		nC = len/nR;

	if (verbose)
		printf("%s - verb(d): %d out(s): %s type(s): float startHigh(d): %d length(%ld) = %d x %d\t wave(s): %s inp(s): %s\n",
				argv[0], verbose, outName, startHigh, len, nC, nR, waveName, inpName);

	if (nC != 0)
	{
		if ( nC * nR > len )
			perror("nC*nR too large, or file too small\n");
		len = nC * nR;
	}

	float 	*floatInput=NULL, *floatOutput=NULL;

	if (verbose)
		printf("%s - verb(d): %d out(s): %s type(s): float startHigh(d): %d length(%ld) = %d x %d\t wave(s): %s inp(s): %s\n",
				argv[0], verbose, outName, startHigh, len, nC, nR, waveName, inpName);

	if (stringNull(outName)) {
		stringCopy(outName, inpName);
	}
	stringCopy(outNameLL, outName);
	stringCopy(outNameLH, outName);
	stringCopy(outNameHL, outName);
	stringCopy(outNameHH, outName);
	strcat(outNameLL, "_LL.fwd");
	strcat(outNameLH, "_LH.fwd");
	strcat(outNameHL, "_HL.fwd");
	strcat(outNameHH, "_HH.fwd");
	
	waveFilterPtr	wfp = wavGetWaveletName  (waveName);
	waveLiftPtr		wlp = wavGetLiftWaveName (waveName);
	assert(wfp);
	assert(wlp);

	if (verbose)
		printf("dwt1d %s: %s --> %s (type float)\n", waveName, inpName, outName);

	FILE* fpInp = fopen(inpName, "r");
	assert (fpInp);
	
	PREPARE_DATA(float);
	fclose(fpInp);
	if (verbose) displayArrayNM_f( floatOutput, nC, nR, stdout );
	
	liftTransform2d_f(floatInput, nC, nR, startHigh, startHigh, false, noRescale, wlp);
	if (verbose) displayArrayNM_f(floatInput, nC, nR, stdout );
				
	float *outLL=NULL, *outLH=NULL, *outHL=NULL, *outHH=NULL;
	long low_nC, low_nR, high_nC, high_nR;
	wavPack_lineBuffers_f ( floatInput, nC, nR, startHigh, startHigh, &outLL, &outHL, &outLH, &outHH, &low_nC, &low_nR, &high_nC, &high_nR);
	if (verbose) 
	{
		printf("\n sizes: low nC nR high nC nR \n %ld %ld %ld %ld \n pointers: %p %p %p %p\n", low_nC, low_nR, high_nC, high_nR, outLL, outHL, outLH, outHH);
		printf("\n\n -------- LL band---------------\n");
		displayArrayNM_f(outLL, low_nC, low_nR, stdout );
		
		printf("\n\n -------- HL band---------------\n");
		displayArrayNM_f(outHL, high_nC, low_nR, stdout );

		printf("\n\n -------- LH band---------------\n");
		displayArrayNM_f(outLH, low_nC, high_nR, stdout );
		
		printf("\n\n -------- HH band---------------\n");
		displayArrayNM_f(outHH, high_nC, high_nR, stdout );
	}
	
	
	FILE* fpOut = fopen(outNameLL, "w");
	assert (fpOut);
	fwrite(outLL, nBytes, low_nC*low_nR, fpOut);
	fclose(fpOut);
	
	fpOut = fopen(outNameHL, "w");
	assert (fpOut);
	fwrite(outHL, nBytes, high_nC*low_nR, fpOut);
	fclose(fpOut);
	
	fpOut = fopen(outNameLH, "w");
	assert (fpOut);
	fwrite(outLH, nBytes, low_nC*high_nR, fpOut);
	fclose(fpOut);

	fpOut = fopen(outNameHH, "w");
	assert (fpOut);
	fwrite(outHH, nBytes, high_nC*high_nR, fpOut);
	fclose(fpOut);

	XYN_DPRINTF(DEBUG_PROCESS,"closed files\n");
	XynFree(floatInput);
	XynFree(floatOutput);
	XynFree(outLL);
	XYN_DPRINTF(DEBUG_PROCESS,"memory freed\n");
	return OK;
}
