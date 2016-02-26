/*
 * xynMem.c
 *
 *  Created on: Mar 11, 2011
 *      Author: mousomer
 */

int xynConvert_i2f ( int* inp, float* out, long len )
{
	assert (inp!=NULL && out!=NULL && len>=0 );
	while (len-->0)
		*out++ = (float) *inp++;

	return OK;
}

int xynConvert_i2c ( int* inp, char* out, long len )
{
	assert (inp!=NULL && out!=NULL && len>=0 );
	while (len-->0)
		*out++ = (char) *inp++;

	return OK;
}

int xynConvert_f2i ( float* inp, int* out, long len )
{
	assert (inp!=NULL && out!=NULL && len>=0 );
	while (len-->0)
		*out++ = (int) *inp++;

	return OK;
}

int xynConvert_f2c ( float* inp, char* out, long len )
{
	assert (inp!=NULL && out!=NULL && len>=0 );
	while (len-->0)
		*out++ = (char) *inp++;

	return OK;
}

int xynConvert_c2i ( char* inp, int* out, long len )
{
	assert (inp!=NULL && out!=NULL && len>=0 );
	while (len-->0)
		*out++ = (int) *inp++;

	return OK;
}

int xynConvert_c2f ( char* inp, float* out, long len )
{
	assert (inp!=NULL && out!=NULL && len>=0 );
	while (len-->0)
		*out++ = (float) *inp++;

	return OK;
}
