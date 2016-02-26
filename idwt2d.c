/*
 * idwt2d.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#include "libxyn.h"

#define USG_STRING "[-v %:] [-out %s:outName] [-type %s:dataType] [-lift %:] [-noScl %:] [-inPlace %:] [-startHigh %:] [-nR %d:nC] %d:nR %s:waveName %s:inpName"
xynString waveName, inpName, outName, dataType;
int verbose=0, lifting=0, inPlace=0, startHigh=0, noRescale=0, nC=0, nR=0;

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
			&verbose, outName, dataType, &lifting, &noRescale, &inPlace, &startHigh, &nC, &nR, waveName, inpName))
		return (FAIL);

	/* transform types:
			0: shuffled,	non lift
			1: inPlace, 	non lift
			2: shuffled, 	lifting implementation
			3: inPlace, 	lifting implementation	*/
	//int transformType = 4*(dataType[0]=='f') + 2*(lifting>0) + (inPlace>0);

	long len = getFileSize(inpName);
	int nBytes = (dataType[0]=='f') ? sizeof(float) : sizeof(int);
	len /= nBytes;
	if (nC == 0)
		nC = len/nR;

	if (verbose)
		printf("%s - verb(d): %d out(s): %s type(s): %s startHigh(d): %d length(%ld) = %d x %d\t wave(s): %s inp(s): %s\n",
				argv[0], verbose, outName, dataType, startHigh, len, nC, nR, waveName, inpName);

	if (nC != 0)
	{
		if ( nC * nR > len )
			perror("nC*nR too large, or file too small\n");
		len = nC * nR;
	}

	int		*intInput=NULL, *intOutput=NULL;
	float 	*floatInput=NULL, *floatOutput=NULL;

	if (verbose)
		printf("%s - verb(d): %d out(s): %s type(s): %s startHigh(d): %d length(%ld) = %d x %d\t wave(s): %s inp(s): %s\n",
				argv[0], verbose, outName, dataType, startHigh, len, nC, nR, waveName, inpName);

	if (stringNull(outName)) {
		stringCopy(outName, inpName);
		strcat(outName, ".fwd");
	}

	if (stringNull(dataType))
		stringCopy(dataType, "d");
	waveFilterPtr	wfp = wavGetWaveletName  (waveName);
	waveLiftPtr		wlp = wavGetLiftWaveName (waveName);
	assert(wfp);
	assert(wlp);

	if (verbose)
		printf("dwt1d %s: %s --> %s (type %c)\n", waveName, inpName, outName, dataType[0]);

	FILE* fpInp = fopen(inpName, "r");
	assert (fpInp);
	FILE* fpOut = fopen(outName, "w");
	assert (fpOut);


	if ( dataType[0] == 'd' ) {
		PREPARE_DATA(int);
		if (verbose) displayArrayNM_i( intOutput, nC, nR, stdout );
		if ( lifting )
			perror ("not implemented. try float\n");
			//liftTransform_i(intInput, len, startHigh, inPlace, wlp);
		else
			perror ("not implemented. try lifting\n");
			//wavTransform_i(intInput, intOutput, len, startHigh, wfp);

		if (verbose) displayArrayNM_i(intOutput, nC, nR, stdout );
		fwrite(intOutput, nBytes, len, fpOut);
	}
	else {
			PREPARE_DATA(float);
			if (verbose) displayArrayNM_f( floatOutput, nC, nR, stdout );
			if ( lifting )
			{
				liftInvTransform2d_f(floatInput, nC, nR, startHigh, startHigh, inPlace, noRescale, wlp);
				if (verbose) displayArrayNM_f(floatInput, nC, nR, stdout );
				fwrite(floatInput, nBytes, len, fpOut);
			}
			else
			{
				perror ("not implemented. try lifting\n");
				//wavInvTransform_f(floatInput, floatOutput, len, startHigh, wfp);
				if (verbose) displayArrayNM_f(floatOutput, nC, nR, stdout );
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
