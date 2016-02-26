/*
 * genRandFile.c
 *
 *  Created on: Jan 31, 2011
 *      Author: mousomer
 */

#include "libxyn.h"
#define USG_STRING "[-lift %:] [-poly %:] %s:waveName"
int lifting = 0, polyPhase = 0;
xynString waveName;

int main(int argc, char* argv[])
{
	if (parseParameters(argc, argv, USG_STRING,
			&lifting, &polyPhase, waveName))
		return (FAIL);

	if (lifting && polyPhase)
		ErrMsg("Either lifting OR polyphase, not both\n.")
	else if (lifting)
	{
		waveLiftPtr lp = wavGetLiftWaveName ( waveName );
		if ( lp == NULL )
			ErrMsg("no wavelet initialized");
		wavLiftPrint ( lp, stdout );
	}
	else if (polyPhase)
	{
		wavePolyPtr lp = wavGetPolyWaveName ( waveName );
		if ( lp == NULL )
			ErrMsg("no wavelet initialized");
		wavPolyPrint ( lp, stdout );
	}
	else
	{
		waveFilterPtr wp = wavGetWaveletName ( waveName );
		if ( wp == NULL )
			ErrMsg("no wavelet initialized");
		wavPrint (wp, stdout);
	}
	
	return OK;
}
