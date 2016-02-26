/*
 * genRandFile.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#include "libxyn.h"
#define BUFFMAX 1024

int main(int argc, char* argv[])
{
	if (argc<3)
		{
			printf ("usage: %s waveName fileName [startHigh] [length]\n", __FUNCTION__);
			return FAIL;
		}


	int nBytes = sizeof(float);
	char inpName[BUFFMAX], outName[BUFFMAX];
	strcpy (inpName, argv[2]);
	strcpy (outName, argv[2]);
	strcat (outName, ".fwd" );

	int startHigh = 0;
	
    long len = getFileSize(inpName) / nBytes;

    if (argc>3)
    {
    	int nVal = 3;
    	if ( strcmp(argv[3], "startHigh") == 0 )
    	{
    		startHigh = 0;
    		if (argc>4)
    			nVal = 4;
    		else nVal = 0;
    	}

    	if (nVal)
    		len = XYN_MIN(atoi(argv[nVal]),len);
    }

    waveletFilterPtr wp = getWaveletName ( argv[1] );
	wavPrint(wp, stdout);

    FILE* fp = fopen ( inpName, "r" );
	assert (fp);
	
	float* data = (float*) malloc ( nBytes * len );
	float* outdata = (float*) malloc ( nBytes * len );
	fread ( data, nBytes, len, fp );
	fclose (fp);

	printf("running function: %s %d %d %s\n", __FUNCTION__, len, startHigh, argv[1]);
	wavTransformFlt( data, outdata, len, startHigh, wp );
	displayFltArray( 	data, len, stdout );
	displayFltArray( outdata, len, stdout );

	fp = fopen ( outName, "w" );
	assert (fp);
	fwrite ( outdata, nBytes, len, fp );
	fclose (fp);
	free(data);
	free(outdata);

	return OK;
}
