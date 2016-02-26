/*
 * genRandFile.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#define MAX_NUM 255
#include "libxyn.h"

#define USG_STRING "[-v %:] [-type %s:dataType] [-min %d:min] [-max %d:max] %d:length %s:fileName"
xynString fName, dataType;
int verbose=0, min=0, max=RAND_MAX, length=10;

int main(int argc, char* argv[])
{
	if (parseParameters(argc, argv, USG_STRING,
				&verbose, &dataType, &min, &max, &length, fName))
			return (FAIL);

/*	char fName[] =  "test123.f";
	char dataType[]="f";*/

	FILE* fp = fopen ( fName, "w" );
	assert (fp);

	int* intData = NULL;
	float* fltData = NULL;

	if (dataType[0]=='d')
	{
		intData = genRandVector_i ( length, max );
		assert (intData);
		XYN_DPRINTF(DEBUG_PROCESS, "data generated: %p\n", intData);
		if (verbose)
			displayArray_i (intData, length, stdout);
		fwrite ( intData, sizeof(int), length, fp );
		XYN_DPRINTF(DEBUG_PROCESS, "data written into %s\n", fName);
		free(intData);
		XYN_DPRINTF(DEBUG_PROCESS, "memory freed\n");
	}
	else
	{
		fltData = genRandVector_f ( length, (float) max );
		assert (fltData);
		XYN_DPRINTF(DEBUG_PROCESS, "data generated: %p\n", fltData);
		if (verbose)
					displayArray_f (fltData, length, stdout);
		fwrite ( fltData, sizeof(int), length, fp );
		XYN_DPRINTF(DEBUG_PROCESS, "data written into %s\n", fName);
		free(fltData);
		XYN_DPRINTF(DEBUG_PROCESS, "memory freed\n");
	}


	fclose (fp);
	XYN_DPRINTF(DEBUG_PROCESS, "file closed.\n");
	return OK;
}
