/*
 * dwt1d.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#include "libxyn.h"

#define USG_STRING "[-v %:] [-out %s:outName] [-noscl %:] [-type %s:dataType] [-lift %:] [-inPlace %:] [-startHigh %:] [-length %d:len] %s:waveName %s:inpName"
xynString waveName, inpName, outName, dataType;
int verbose=0, lifting=0, inPlace=0, startHigh=0, length=0, noRescale=0;

#define PREPARE_DATA(_t_) {		\
		nBytes = sizeof(_t_);						\
		len = (length > 0) ? length : len/nBytes;	\
		XYN_DPRINTF(DEBUG_PROCESS,"\npreparing data. Length %d/=%d\n", len, nBytes);	\
		_t_##Input  = (_t_*) malloc(nBytes * len);	\
		_t_##Output = (_t_*) malloc(nBytes * len);	\
		assert(_t_##Input && _t_##Output );			\
		fread(_t_##Input, nBytes, len, fpInp);		\
		memcpy(_t_##Output, _t_##Input, len*nBytes );		\
		if (verbose) displayArray_##_t_(_t_##Output, len, stdout );	 }


int main(int argc, char* argv[])
{
	if (parseParameters(argc, argv, USG_STRING,
			&verbose, outName, &noRescale, dataType, &lifting, &inPlace, &startHigh, &length, waveName, inpName))
		return (FAIL);

	/* transform types:
			0: shuffled,	non lift
			1: inPlace, 	non lift
			2: shuffled, 	lifting implementation
			3: inPlace, 	lifting implementation	*/
	//int transformType = 4*(dataType[0]=='f') + 2*(lifting>0) + (inPlace>0);

	if (verbose)
		printf("%s - verb(d): %d out(s): %s type(s): %s startHigh(d): %d length(d): %d wave(s): %s inp(s): %s\n",
				argv[0], verbose, outName, dataType, startHigh, length, waveName, inpName);

	if (stringNull(outName)) {
		stringCopy(outName, inpName);
		strcat(outName, ".fwd");
	}

	if (stringNull(dataType))
		stringCopy(dataType, "d");
	waveFilterPtr	wfp = wavGetWaveletName  (waveName);
	waveLiftPtr		wlp = wavGetLiftWaveName (waveName);
	assert(wfp); 	assert(wlp);

	if (verbose)
		printf("dwt1d %s: %s --> %s (type %c)\n", waveName, inpName, outName, dataType[0]);

	FILE* fpInp = fopen(inpName, "r");
	assert (fpInp);
	FILE* fpOut = fopen(outName, "w");
	assert (fpOut);

	long len = getFileSize(inpName);
	int nBytes = sizeof(int);
	int		*intInput=NULL, *intOutput=NULL;
	float 	*floatInput=NULL, *floatOutput=NULL;

	if ( dataType[0] == 'd' ) {
		PREPARE_DATA(int);

		if ( !noRescale )
			upScaleVect_i( intInput, len );

		if ( lifting )
		{
			liftTransform_i(intInput, len, startHigh, inPlace, noRescale, wlp);
		//	rescaleVect_i(intInput, len);

			if (verbose) displayArray_i(intInput, len, stdout );
			fwrite(intInput, nBytes, len, fpOut);
		}
		else
		{
			wavTransform_i(intInput, intOutput, len, startHigh, wfp);
		//	rescaleVect_i(intInput, len);

			if (verbose) displayArray_i(intOutput, len, stdout );
			fwrite(intOutput, nBytes, len, fpOut);
		}
	}
	else {
			PREPARE_DATA(float);

			if ( lifting )
			{
				liftTransform_f(floatInput, len, startHigh, inPlace, noRescale, wlp);
				if (verbose) displayArray_f(floatInput, len, stdout );
				fwrite(floatInput, nBytes, len, fpOut);
			}
			else
			{
				wavTransform_f(floatInput, floatOutput, len, startHigh, inPlace, wfp);
				if (verbose) displayArray_f(floatOutput, len, stdout );
				fwrite(floatOutput, nBytes, len, fpOut);
			}
		}

	
	fclose(fpInp);fclose(fpOut);
	XYN_DPRINTF(DEBUG_PROCESS,"closed files\n");
	XynFree(intInput);	XynFree(floatInput);
	XynFree(intOutput);	XynFree(floatOutput);
	XYN_DPRINTF(DEBUG_PROCESS,"memory freed\n");
	return OK;
}
