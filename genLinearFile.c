/*
 * genRandFile.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#include "libxyn.h"
#define USG_STRING "[-v %:] [-type %s:dataType] [-power %d:power] [-min %d:min] %d:max %d:len %s:fileName"
xynString fName, dataType;
int verbose=0, power=1, min=0, max=0, len=0;

int main(int argc, char* argv[])
{
	if (parseParameters(argc, argv, USG_STRING, &verbose,
			dataType, &power, &min, &max, &len, fName))
		return (FAIL);

	if (verbose)
		printf("%s - generating array into file %s. data of type %c, from %d to %d length %d. ",
				__FILE__, fName, dataType[0], min, max, len);

	FILE* fp = fopen ( fName, "w" );
	assert (fp);

	int nBytes = sizeof(float);
	float* data = (float*) malloc (len * nBytes );

	switch (power)
	{
		case 0:
			if (verbose)
				printf(" Constant array.\n");
			constVect_f(data, len, min, max);
			break;
		case 1:
			if (verbose)
				printf(" Linear array.\n");
			linearVect_f(data, len, min, max);
			break;
		default:
			if (verbose)
				printf(" Polynomial array of order %d.\n", power);
			powerVect_f(data, len, min, max, power);
	}
	assert (data);

	if (verbose)
		displayArray_f(data, len, stdout);

	if (dataType[0]=='d')
	{
		nBytes = sizeof(int);
		int* intData = (int*) malloc (len * nBytes);
		convert_ftoi ( intData, data, len );
		fwrite ( intData, nBytes, len, fp );
		free ( intData );
	}
	else if (dataType[0]=='g')
	{
		nBytes = sizeof(double);
		double* doubleData = (double*) malloc (len * nBytes);
		convert_ftog( doubleData, data, len );
		fwrite ( doubleData, nBytes, len, fp );
		free ( doubleData );
	}
	else
		fwrite ( data, nBytes, len, fp );

	fclose (fp);
	free(data);

	return 0;
}
