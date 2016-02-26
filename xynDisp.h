/*
 * xynDisp.h
 *
 *  Created on: Feb 6, 2011
 *      Author: mousomer
 */

#ifndef XYNDISP_H_
#define XYNDISP_H_

inline void XYN_DPRINTF(int level, char *format, ...)
{
#    ifdef DEBUG
	va_list args;
	va_start(args, format);
	if(debugPrint & level) {
		vfprintf(stderr, format, args);
	}
	va_end(args);
#    endif /* _DEBUG */
}

//void DPRINTF(int level, char *format, ...);

void disp_i ( int* num);
void disp_c ( char* num);
void disp_f ( float* num);
void disp_e ( double* num);
int displayArray ( char* arr, int nBytes, int len, void(*dispNumber)(void*) );
int displayArrayFormat ( void* arr, int nBytes, int len, char type, int digits, int precision, int lineLen, FILE* fp );

int displayArray_i ( int* arr, int len, FILE* fp );
#define displayArray_int displayArray_i
int displayArray_f ( float* arr, int len, FILE* fp );
int displayArray_d ( double* arr, int len, FILE* fp );

#define displayArray_float displayArray_f
int displayArrayJump_f ( float* arr, int jmp, int len, FILE* fp );
int displayArrayJump_d ( double* arr, int jmp, int len, FILE* fp );

int displayArrayNM_i ( int* arr, int nC, int nR, FILE* fp );
#define displayArrayNM_int displayArray_i
int displayArrayNM_f ( float* arr, int nC, int nR, FILE* fp );
int displayArrayNM_d ( double* arr, int nC, int nR, FILE* fp );
#define displayArrayNM_float displayArray_f

#endif /* XYNDISP_H_ */
