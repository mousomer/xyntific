/*
 * lift_dwt1d.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#include "libxyn.h"
/*extern struct waveletFilter wavCDF53;
extern struct waveletFilter wavCDF97;
extern struct waveletFilter wavHaar;*/

#define USG_STRING "[-V %:] [-out %s:outName] [-type %s:dataType] [-startHigh %:] [-length %d:len] %s:waveName %s:inpName"
xynString waveName, inpName, outName, dataType;
int verbose=0, startHigh=0, length=0;

int main(int argc, char* argv[])
{
	if (parseParameters(argc, argv, USG_STRING, &verbose, outName,
						dataType, &startHigh, &length, waveName, inpName))
		return (FAIL);

	if (verbose)
		printf("%s - verb(d): %d out(s): %s type(s): %s startHigh(d): %d length(d): %d wave(s): %s inp(s): %s\n",
				argv[0], verbose, outName, dataType, startHigh, length, waveName, inpName);
	if (stringNull(outName)) {
		stringCopy(outName, inpName);
		strcat(outName, ".fwd");
	}

	if (stringNull(dataType))
		stringCopy(dataType, "d");
	waveFilterPtr wp = getWaveletName(waveName);
	assert(wp);

	if (verbose)
	{
		printf("dwt1d %s: %s --> %s (type %c)\n", waveName, inpName, outName, dataType[0]);
		wavPrint(wp, stdout);
	}

	FILE* fpInp = fopen(inpName, "r");
	assert (fpInp);
	FILE* fpOut = fopen(outName, "w");
	assert (fpOut);

	long len = getFileSize(inpName);
	int nBytes = sizeof(int);
	int		*intInput=NULL, *intOutput=NULL;
	float 	*fltInput=NULL, *fltOutput=NULL;
	long 	*lngInput=NULL, *lngOutput=NULL;
	double	*dblInput=NULL, *dblOutput=NULL;

	switch (dataType[0]) {
	case 'd':
		nBytes = sizeof(int);
		len = (length > 0) ? length : len/nBytes;
		intInput  = (int*) malloc(nBytes * len);
		intOutput = (int*) malloc(nBytes * len);
		assert(intInput && intOutput );
		fread(intInput, nBytes, len, fpInp);
		if (verbose) displayIntArray(intInput, len, stdout );
		wavTransform_d(intInput, intOutput, len, startHigh, wp);
		if (verbose) displayIntArray(intOutput, len, stdout );
		fwrite(intOutput, nBytes, len, fpOut);
		break;

	case 'f':
		nBytes = sizeof(float);
		len = (length > 0) ? length : len/nBytes;
		printf("%dx%d\n", nBytes, len);
		fltInput  = (float*) malloc(nBytes * len);
		fltOutput = (float*) malloc(nBytes * len);
		assert(fltInput && fltOutput );
		fread(fltInput, nBytes, len, fpInp);
		if (verbose) displayFltArray(fltInput, len, stdout );
		wavTransform_f(fltInput, fltOutput, len, startHigh, wp);
		if (verbose) displayFltArray(fltOutput, len, stdout );
		fwrite(fltOutput, nBytes, len, fpOut);
		break;

	default:
		printf ("type not supported\n");
	}
	
	fclose(fpInp);fclose(fpOut);
	XynFree(intInput);	XynFree(fltInput);
	XynFree(dblInput);	XynFree(lngInput);
	XynFree(intOutput);	XynFree(fltOutput);
	XynFree(dblOutput);	XynFree(lngOutput);

	return OK;
}
