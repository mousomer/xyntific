/*
 * bkup.c
 *
 *  Created on: Mar 14, 2011
 *      Author: mousomer
 */


function fd()
{
		float* tmp = (float*) malloc (len * sizeof(float));
		memcpy (tmp, data, len*sizeof(float));
		float* tmpCpy = tmp;
		long lowLen = len/2 + (len%2)*sclStart;
		//long highLen = len - lowLen;
		float* outLH[2] = {data, data + lowLen};

		int lh = wavStart;
		for (long i=0; i<len; i++)
		{
			*(outLH[lh])++ = *tmpCpy++;
			lh = 1-lh;
		}
		free(tmp);
}

int wavTransformInPlace_f( float* dataInp, float* dataOut, long len, bool oddity, waveFilterPtr filt )
{
	assert(dataInp!=NULL); assert(dataOut!=NULL); assert(filt!=NULL);
	assert(len>0);

	float  *cpyOut[2] = { dataOut+oddity, dataOut+1-oddity }, 	*cpyInp = dataInp;

	int boundary_down = MAX( filt->pos_anal[0]-oddity, filt->pos_anal[1]-(1-oddity));
	int boundary_up = len - boundary_down;
	XYN_DPRINTF(DEBUG_HEADER, "%s %s: len %d, oddity %d\n", __FUNCTION__, filt->name, len, oddity);
	int i=0;
	for (; i<boundary_down; i++, oddity=1-oddity)
	{
		*(cpyOut[oddity]) = dotProd_f ( filt->f_anal_folded[oddity][i], cpyInp, (filt->len_anal[oddity])/2 +1 );
		(cpyOut[oddity])+=2;
	}

	int fold=i+1;
	for (; i<boundary_up; i++, oddity=1-oddity)
	{
		*(cpyOut[oddity]) = dotProd_f ( filt->f_anal[oddity], cpyInp++, filt->len_anal[oddity] );
		(cpyOut[oddity])+=2;
	}

	for (; i<len; i++, fold++, oddity=1-oddity)
	{
		*(cpyOut[oddity]) = dotProd_f ( filt->f_anal_folded[oddity][fold], cpyInp, (filt->len_anal[oddity])/2 +1 );
		(cpyOut[oddity])+=2;
	}

	return OK;
}

static int wavInvTransformHiLow_f( float* dataOut, long len, float* inpLow, long lenLow,
							float* inpHigh, long lenHigh, bool oddity, waveFilterPtr filt )
{
	assert(dataOut!=NULL); assert(inpLow!=NULL); assert(inpHigh!=NULL); assert(filt!=NULL);
	assert(lenLow+lenHigh==len);
	float  *cpyInp[2] = { inpLow, inpHigh }, 	*cpyOut = dataOut;
	int boundary_down = MAX( filt->pos_synt[0]-oddity, filt->pos_synt[1]-(1-oddity));
	int boundary_up = len - boundary_down;

	XYN_DPRINTF(DEBUG_HEADER, "%s %s: len %d=%d+%d, (boundaries %d->%d) oddity %d\n",
			__FUNCTION__, filt->name, len, lenLow, lenHigh, boundary_down, boundary_up, oddity);

	int i=0;
	for (; i<boundary_down; i++, oddity=1-oddity)
	{
		*cpyOut  = dotProdJump_f ( filt->f_synt_folded[oddity][i], 2, cpyInp[oddity], 1, (filt->len_synt[oddity])/2 +1 );
		oddity=1-oddity;
		*cpyOut += dotProdJump_f ( filt->f_synt_folded[oddity][i], 2, cpyInp[oddity], 1, (filt->len_synt[oddity])/2 +1 );
		cpyOut++;
	}

	displayArray_f(filt->f_synt[oddity], filt->len_synt[oddity], stdout);
	displayArray_f(filt->f_synt[1-oddity], filt->len_synt[1-oddity], stdout);
	int fold=i+1;
	for (; i<boundary_up; i++)
	{
		*cpyOut  = dotProdJump_f ( filt->f_synt[oddity], 2, cpyInp[oddity], 1, (filt->len_synt[oddity])/2 +1 );
		oddity = 1-oddity;
		*cpyOut += dotProdJump_f ( filt->f_synt[oddity], 2, cpyInp[oddity], 1, (filt->len_synt[oddity])/2 +1 );
		cpyOut++;
	}

	for (; i<len; i++, fold++)
	{
		*cpyOut  = dotProdJump_f ( filt->f_synt_folded[oddity][fold], 2, cpyInp[oddity], 1, (filt->len_synt[oddity])/2 +1 );
		oddity=1-oddity;
		*cpyOut += dotProdJump_f ( filt->f_synt_folded[oddity][fold], 2, cpyInp[oddity], 1, (filt->len_synt[oddity])/2 +1 );
		cpyOut++;
	}

	return OK;
}


/* floating point transform*/
int wavInvTransformInPlace_f( float* dataInp, float* dataOut, long len, bool oddity, waveFilterPtr filt )
{
	assert(dataInp!=NULL); assert(dataOut!=NULL); assert(filt!=NULL);
	assert(len>0);

	float  *cpyOut[2] = { dataOut+oddity, dataOut+1-oddity }, 	*cpyInp = dataInp;

	int boundary_down = MAX( filt->pos_synt[0]-oddity, filt->pos_synt[1]-(1-oddity));
	int boundary_up = len - boundary_down;

	int i=0;
	for (; i<boundary_down; i++, oddity=1-oddity)
	{
		*(cpyOut[oddity]) = dotProd_f ( filt->f_synt_folded[oddity][i], cpyInp, (filt->len_synt[oddity])/2 +1 );
		(cpyOut[oddity])+=2;
	}

	int fold=i+1;
	for (; i<boundary_up; i++, oddity=1-oddity)
	{
		*(cpyOut[oddity]) = dotProd_f ( filt->f_synt[oddity], cpyInp++, filt->len_synt[oddity] );
		(cpyOut[oddity])+=2;
	}

	for (; i<len; i++, fold++, oddity=1-oddity)
	{
		*(cpyOut[oddity]) = dotProd_f ( filt->f_synt_folded[oddity][fold], cpyInp, (filt->len_synt[oddity])/2 +1 );
		(cpyOut[oddity])+=2;
	}

	return OK;
}


int wavTrasformHi-Low ( float* dataOut, long len, float* inpLow, long lenLow )
{
	/* oddity boundary (lowpass if oddity=0, else highpass */
	int j=oddity;
	for (i=0; i<boundary_down[j]; i+=2)		// i from 0 to boundary - this is first pass type
	{

		*(cpyOut[j]) = dotProd_f ( filt->f_anal_folded[j][fold[j]], cpyInp[j], filt->len_anal[j] );
		XYN_DPRINTF(DEBUG_DETAIL, "i=%6ld; fold %4d=%7.2f(oddity %2d)\t", i, fold[j], *cpyOut[j], j);
		displayArray_f(filt->f_anal_folded[j][fold[j]], filt->len_anal[j], stdout);

		cpyOut[j] ++ ;								// output pointer advance
		fold[j]+=2;
	}
	XYN_DPRINTF(DEBUG_PROCESS, "%s lower boundary done for %d. i=%d, folds (%d %d)\n", __FUNCTION__, j, i, fold[0], fold[1]);
	/* 1-oddity boundary */
	j=1-oddity;
	for (i=1; i<boundary_down[j]; i+=2)		// i from 1 - this is second pass type
	{
		*(cpyOut[j]) = dotProd_f ( filt->f_anal_folded[j][fold[j]], cpyInp[j], filt->len_anal[j] );
		XYN_DPRINTF(DEBUG_DETAIL, "i=%6ld; fold %4d=%7.2f(oddity %2d)\t", i, fold[j], *cpyOut[j], j);
		displayArray_f(filt->f_anal_folded[j][fold[j]], filt->len_anal[j], stdout);

		cpyOut[j] ++ ;								// output pointer advance
		fold[j]+=2;
	}
	XYN_DPRINTF(DEBUG_PROCESS, "%s lower boundary done for %d. i=%d, folds (%d %d)\n", __FUNCTION__, j, i, fold[0], fold[1]);

	/* small boundary fill */
	int diffBound = abs ( filt->pos_anal[oddity] - filt->pos_anal[1-oddity] );		// small boundary
	bool j = (filt->pos_anal[oddity] < filt->pos_anal[1-oddity]) ? oddity : 1-oddity;
	for (i=(j!=oddity); i<boundary_down[1-j]; i+=2 )
	{
		*(cpyOut[j]) = dotProd_f ( filt->f_anal[j], cpyInp[j], filt->len_anal[j] );
		XYN_DPRINTF(DEBUG_DETAIL, "- %3ld=%6.2f(%d)\n", i, *cpyOut[j], j); if (i%8==7) NL;
		cpyOut[j] += 1;
		cpyInp[j] += 2;
	}



	if ( fold[oddity] - filt->pos_anal[oddity] == 1)		// input pointer - advance by 1 if exceeded boundary by 1
		cpyInp[oddity] ++;
	else if ( fold[oddity] == filt->pos_anal[oddity] )		// folding pointer - advance by 2 into upper boundary
		fold[oddity] +=2;
	if ( fold[1-oddity] - filt->pos_anal[1-oddity] == 1)	// input pointer - advance by 1 if exceeded boundary by 1
		cpyInp[1-oddity] ++;
	else if ( fold[1-oddity] == filt->pos_anal[1-oddity] )		// folding pointer - advance by 2 into upper boundary
		fold[1-oddity] +=2;

	/* main loop - run between boundaries, input incremented by 2 for each pass-type */
	for (; i<boundary_up; i++, oddity=1-oddity)
	{
		*(cpyOut[oddity]) = dotProd_f ( filt->f_anal[oddity], cpyInp[oddity], filt->len_anal[oddity] );
		//printf("- %3ld=%6.2f(%d)\t", i, *cpyOut[oddity], oddity); if (i%8==7) NL;
		cpyOut[oddity] += 1;
		cpyInp[oddity] += 2;
	}

	//fold[0]++; fold[1]++;
	for (; i<len; i++, oddity=1-oddity)
	{

		if ( len-i > filt->pos_anal[oddity] )
		{
			*(cpyOut[oddity])++ = dotProd_f ( filt->f_anal[oddity], cpyInp[oddity], (filt->len_anal[oddity]) );
			XYN_DPRINTF(DEBUG_DETAIL, "i=%6ld; not folded =%7.2f(oddity %2d)\n", i, *cpyOut[oddity], oddity);
			//displayArray_f(filt->f_anal_folded[oddity][fold[oddity]], filt->len_anal[oddity], stdout);
			cpyInp[oddity] += 2;
		}
		else
		{
			*(cpyOut[oddity])++ = dotProd_f ( filt->f_anal_folded[oddity][fold[oddity]], cpyInp[oddity], filt->len_anal[oddity] );
			XYN_DPRINTF(DEBUG_DETAIL, "i=%6ld; fold %4d=%7.2f(oddity %2d)\t", i, fold[oddity], *cpyOut[oddity], oddity);
			displayArray_f(filt->f_anal_folded[oddity][fold[oddity]], filt->len_anal[oddity], stdout);
			fold[oddity]+=2;
		}
	}

}

