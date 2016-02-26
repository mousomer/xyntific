/*
 * genRandFile.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */
#include "libxyn.h"


#define USG_STRING "[-prec %d:precision] [-digits %d:digits] [-line %d:lineLen] [-type %: %s:dataType(d|c|l|f|g)] [-length %: %d:len] %s:fName"
xynString fName, dataType;
int precision=2, digits=6, lengthSpecified=0, typeSpecified=0, length=0, lineLen=0;

int main(int argc, char* argv[])
{
	  if (parseParameters(argc, argv, USG_STRING,
	                         &precision, &digits, &lineLen, &typeSpecified, dataType,
	                         &lengthSpecified, &length, fName))  return(FAIL);

	int nBytes = 0;
	if ( !typeSpecified ) dataType[0]='d';
	
	switch (dataType[0])
	{
	case 'f': 
		nBytes = sizeof(float);
		printf("printing float data from %s\n", fName);
		break;
	case 'g':
		nBytes = sizeof(double);
		printf("printing double data from %s\n", fName);
		break;
	case 'e':
		nBytes = sizeof(double);
		printf("printing double data from %s\n", fName);
		break;
	case 'c':
		nBytes = sizeof(char);
		printf("printing char data from %s\n", fName);
		break;
	default:
		nBytes = sizeof(int);
		printf("printing integer data from %s\n", fName);
		break;
	}

    long len = getFileSize(fName) / nBytes;
    if ( lengthSpecified )
    	len = XYN_MIN(len, length);

	FILE* fp = fopen ( fName, "r" );
	assert (fp);
	
	char* data = malloc ( nBytes * len );
	printf( "pointer %p   requested %ld=%d x %ld bytes \n", data, nBytes*len, nBytes, len);
	fread ( data, nBytes, len, fp );

	//displancyArrayFmt ( data, nBytes, len, format, stdout );
	//void* dataCpy = data;
	displayArrayFormat ( data, nBytes, len, dataType[0], digits, precision, lineLen, stdout );
	printf( "data is %ld long (%d bytes each)  \n", len, nBytes);

	fclose (fp);
	free(data);

	return 0;
}
