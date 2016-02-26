/*
 * calcMoment.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#include "libxyn.h"

#define USG_STRING "[-v %:] [-type %s:dataType] [-central %:] [-divImp %:] [-length %d:len] %d:moment %s:inpName"
xynString inpName, dataType;
int verbose=0, central=0, moment=0, divImp=0, length=0;

int main(int argc, char* argv[])
{
	if (parseParameters(argc, argv, USG_STRING,
			&verbose, dataType, &central, &divImp, &length, &moment, inpName))
		return (FAIL);

	/* calc moment types:
			0: integer,	non central, 	not normed
			1: integer,	central, 		not normed
			2: integer,	not central, 	normed
			3: integer,	central, 		normed

			5: float,	non central, 	not normed
			6: float,	central, 		not normed
			7: float,	not central, 	normed
			8: float,	central, 		normed

	int calcType =  (central) ? 1 : 0;
		calcType += (normed)  ? 2 : 0;
		//calcType += (dataType[0] == 'f') ? 5 : 0;*/

	if (stringNull(dataType))
		stringCopy(dataType, "d");

	if (verbose)
		printf("%s - verb: %d, dataType: %s, central: %d, divImp: %d, moment: %d, length: %d, inp(s): %s\n",
				argv[0], verbose, dataType, central, divImp, moment, length, inpName);

	FILE* fpInp = fopen(inpName, "r");
	assert (fpInp);

	long len = getFileSize(inpName);
	XYN_DPRINTF(DEBUG_DETAIL,"file size %ld\n", len);
	int nBytes = 1;
	int		*intInput=NULL;
	float 	*fltInput=NULL;

	double out = 0;

	switch (dataType[0]) {
	case 'd':
		nBytes = sizeof(int);
		len = ((length>0) && (length<len/nBytes)) ? length : len/nBytes;
		intInput  = (int*) malloc(nBytes * len);
		assert( intInput );
		fread(intInput, nBytes, len, fpInp);
		if (verbose) displayArray_i(intInput, len, stdout );

		if ( divImp )
			out = momentDiv_i( intInput, len, moment, central );
		else
			out = moment_i	( intInput, len, moment, central );
		break;

	case 'f':
		nBytes = sizeof(float);
		len = ((length>0) && (length<len/nBytes)) ? length : len/nBytes;
		fltInput  = (float*) malloc(nBytes * len);
		assert(fltInput);
		fread(fltInput, nBytes, len, fpInp);
		if (verbose) displayArray_f(fltInput, len, stdout );
		XYN_DPRINTF(DEBUG_DETAIL,"data length %ld /= %d\n", len, nBytes );

		if ( divImp )
			out = momentDiv_f	( fltInput, len, moment, central );
		else
			out = moment_f	( fltInput, len, moment, central );
		break;

	default:
		printf ("type not supported\n");
	}

	printf("moment %d_th is %e\n", moment, out);
	fclose(fpInp);
	XYN_DPRINTF(DEBUG_DETAIL,"closed files\n");
	XynFree(intInput);	XynFree(fltInput);
	XYN_DPRINTF(DEBUG_DETAIL,"memory freed\n");
	return OK;
}

