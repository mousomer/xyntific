/*
 * txt2bin.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */
#include "libxyn.h"

#define MAX_BUFF 1024
#define USG_STRING "[-type %: %s:dataType(d|c|l|f|e)] [-length %: %d:len] [-outName %: %s:outName] %s:fileName"
xynString inpName, outName, dataType;
int length=MAX_BUFF, typeSpecified=0, lengthSpecified=0, outNameSpecified=0;

int main(int argc, char* argv[])
{
	  if (parseParameters(argc, argv, USG_STRING,
	                         &typeSpecified, dataType, &lengthSpecified, &length, &outNameSpecified, outName, inpName))
		  return(FAIL);

	int nBytes=0, len=0;
	int *i_data=NULL;
	float *f_data=NULL;
	if ( !typeSpecified ) dataType[0]='d';
	if ( !outNameSpecified )
	{
		stringCopy(outName, inpName);
		
	}
	FILE* fp = fopen ( inpName, "r" );
	assert (fp);

	// read file length
	if ( getIntFromFile ( fp, 1, &len ) != OK )
		ErrMsg ( "trouble reading segment size from file.\n");
	length = (length<len) ? length : len;

	if (dataType[0]=='f')
	{
		strcat(outName, ".f");
		nBytes = sizeof(float);
		printf("extracting float data from %s\n", inpName);
		f_data = malloc ( nBytes * len );

		for (int i=0; i<length; i++)
		{
			if ( getFloatFromFile ( fp, 1, &(f_data[i]) ) != OK )
						ErrMsg ( "trouble reading subpixel rsolution from file.\n" );
			XYN_DPRINTF(DEBUG_PROCESS," = %f\n", f_data[i] );
		}

		fclose(fp);
		fp = fopen (outName, "w");
		fwrite ( f_data, nBytes, length, fp );
		free(f_data);
	}
	else
	{
		strcat(outName, ".i");
		nBytes = sizeof(int);
		printf("extracting integer data from %s\n", inpName);
		i_data = malloc ( nBytes * len );

		for (int i=0; i<length; i++)
		{
			if ( getIntFromFile ( fp, 1, &(i_data[i]) ) != OK )
						ErrMsg ( "trouble reading subpixel rsolution from file.\n" );
			XYN_DPRINTF(DEBUG_PROCESS," = %ld\n", i_data[i]);
		}

		fclose(fp);
		fp = fopen (outName, "w");
		fwrite ( i_data, nBytes, length, fp );
		free(i_data);
	}

	
	fclose (fp);

	return OK;
}
