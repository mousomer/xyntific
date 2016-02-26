/*
 * xynTmpWave.c
 *
 *  Created on: Mar 15, 2011
 *      Author: mousomer
 *
 *      temporal transform on wavelet
 */
#include "libxyn.h"
static int FixWaveShift = 10;
static int FixWaveScale = 1000;



tempBuffPtr pBuffer;
int initMacroBlockBuffer ( segmInfPtr sp );
int* getMacroBlockFrame ( XYNCh)
static int updateBuffer ( FILE* inpFp, segmInfPtr sp, int frame );



static int updateBuffer ( FILE* inpFp, segmInfPtr sp, int frame )
{
	/* sanity check */
	assert (inpFp!=NULL); assert (sp!=NULL); assert (frame>=0 && frame<sp->size[2]);
	if ( frame)

}






/* end of file */
