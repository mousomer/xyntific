/*
 * xynDisp.c
 *
 *  Created on: Feb 6, 2011
 *      Author: mousomer
 */

#include "libxyn.h"
#define PRNT_LINE_LEN 10

#define PRINT_FMT_ARRAY(fmt,digits,y)	\
	for ( int i=0; i<len; i++ )	{		\
		fprintf(fp, fmt, digits, y );	\
		if ( (i+1)%lineLen == 0 ) NL;}

#define PRINT_FMT_PREC_ARRAY(fmt,digits, precision,y)	\
	for ( int i=0; i<len; i++ )	{						\
		fprintf(fp, fmt, precision, digits, y );		\
		if ( (i+1)%lineLen == 0 ) NL;	}


/* debugger printing information */
/*int debugPrint = DEBUG_FUNCTION + DEBUG_HEADER + DEBUG_TIME +DEBUG_COEFF + DEBUG_PROCESS;*/
int debugPrint = DEBUG_FUNCTION + DEBUG_HEADER + DEBUG_TIME +DEBUG_COEFF + DEBUG_PROCESS + DEBUG_DETAIL;
extern inline void XYN_DPRINTF(int level, char *format, ...);


void disp_i( int* num)
{
	printf(" %6d ", *num );
}
void disp_c( char* num)
{
	printf(" %4c ", *num );
}
void disp_f( float* num)
{
	printf(" %f ", *num );
}
void disp_e( double* num)
{
	printf(" %lf ", *num );
}

int displayArray ( char* arr, int nBytes, int len, void(*dispNumber)(void*) )
{
	assert ( arr!=NULL && nBytes>0 && len>=0 );

	printf ("\nDisplaying array... \n");
	while ( len-- > 0 )
	{
		if ( (len+1)%PRNT_LINE_LEN == 0 ) NL;
		dispNumber ( arr );
		arr += nBytes;
	}

	NL;
	return OK;
}

int displayArrayFormat ( void* arr, int nBytes, int len, char type, int digits, int precision, int lineLen, FILE* fp )

{
	assert ( arr!=NULL && nBytes>0 && len>=0 );
	XYN_DPRINTF(DEBUG_HEADER,"\nDisplaying array... %p, %dX%d, type %c, %d.%d at %d \n", arr, nBytes, len, type, digits, precision, lineLen);

	int* int_arr;
	char* char_arr;
	float* float_arr;
	long* long_arr;
	double* double_arr;

	if ( !lineLen ) lineLen = PRNT_LINE_LEN;
	switch (type)
	{
	case 'd':
		int_arr = (int*) (arr);
		PRINT_FMT_ARRAY(" %*d ", digits, *int_arr++);
		break;
	case 'l':
		long_arr = (long*) (arr);
		PRINT_FMT_ARRAY(" %*ld ", digits, *long_arr++);
		break;
	case 'c':
		char_arr = (char*) (arr);
		PRINT_FMT_ARRAY(" %*c ", digits, *char_arr++);
		break;
	case 'f':
		float_arr = (float*) (arr);
		PRINT_FMT_PREC_ARRAY(" %*.*f ", precision, digits, *float_arr++);
		break;
	case 'g':
		double_arr = (double*) (arr);
		PRINT_FMT_PREC_ARRAY(" %*.*lg ", precision, digits, *double_arr++);
		break;
	default:
		perror ( "data type unknown ");
	}
	NL;
	return OK;
}

int displayArray_i ( int* arr, int len, FILE* fp )
{
	assert ( arr!=NULL && fp && len>=0 );

	while ( len-- > 0 )
	{
		if ( (len+1)%PRNT_LINE_LEN == 0 ) NL;
		fprintf(fp, " %5d ", *arr++);
	}
	NL;
	return OK;
}

int displayArray_d ( double* arr, int len, FILE* fp )
{
	assert ( arr!=NULL && fp && len>=0 );

	while ( len-- > 0 )
	{
		if ( (len+1)%PRNT_LINE_LEN == 0 ) NL;
		fprintf(fp, " %12.4f ", *arr++);
	}
	NL;
	return OK;
}

int displayArray_f ( float* arr, int len, FILE* fp )
{
	assert ( arr!=NULL && fp && len>=0 );

	while ( len-- > 0 )
	{
		if ( (len+1)%PRNT_LINE_LEN == 0 ) NL;
		fprintf(fp, " %12.4f ", *arr++);
	}
	NL;
	return OK;
}

int displayArrayJump_f ( float* arr, int jmp, int len, FILE* fp )
{
	assert ( arr!=NULL && fp && len>=0 );

	while ( len-- > 0 )
	{
		if ( (len+1)%PRNT_LINE_LEN == 0 ) NL;
		fprintf(fp, " %12.4f ", *arr);
		arr += jmp;
	}
	NL;
	return OK;
}

int displayArrayJump_d ( double* arr, int jmp, int len, FILE* fp )
{
	assert ( arr!=NULL && fp && len>=0 );

	while ( len-- > 0 )
	{
		if ( (len+1)%PRNT_LINE_LEN == 0 ) NL;
		fprintf(fp, " %12.4f ", *arr);
		arr += jmp;
	}
	NL;
	return OK;
}

int displayArrayNM_i ( int* arr, int nC, int nR, FILE* fp )
{
	assert ( arr!=NULL && fp && nC>=0 && nR>=0 );

	while ( nR-- > 0 )
	{
		displayArray_i ( arr, nC, fp );
		arr+=nC;
	}
	NL;
	return OK;
}

int displayArrayNM_f ( float* arr, int nC, int nR, FILE* fp )
{
	assert ( arr!=NULL && fp && nC>=0 && nR>=0 );

	while ( nR-- > 0 )
	{
		displayArray_f ( arr, nC, fp );
		arr+=nC;
	}
	NL;
	return OK;
}

int displayArrayNM_d ( double* arr, int nC, int nR, FILE* fp )
{
	assert ( arr!=NULL && fp && nC>=0 && nR>=0 );

	while ( nR-- > 0 )
	{
		displayArray_d ( arr, nC, fp );
		arr+=nC;
	}
	NL;
	return OK;
}
