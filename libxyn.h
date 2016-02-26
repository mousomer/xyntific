/*
 * libxyn.h
 *
 *  Created on: 07/01/2011
 *      Author: lilah215
 *
 *
 *
 *      todos:
 *
 *      	1. get non-inplace checked
 *      	2. get 2d + 2d cropped versions of dwt -> dwt1d with jumps
 *      	3. do dwt1d on frames
 *
 *
 */

#ifndef LIBXYN_H_
#define LIBXYN_H_

//#define DEBUG
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdbool.h>


#define DEBUG_FUNCTION 	(1<< 0)
#define DEBUG_HEADER	(1<< 1)
#define DEBUG_PROCESS	(1<< 2)
#define DEBUG_COEFF		(1<< 3)
#define DEBUG_DETAIL 	(1<< 4)
#define DEBUG_TIME 		(1<< 5)

extern int debugPrint; // = DEBUG_FUNCTION+DEBUG_HEADER+DEBUG_DETAIL+DEBUG_TIME;


/* type sizes */
static const int LNG_BIT = CHAR_BIT * sizeof(long);
static const int INT_BIT = CHAR_BIT * sizeof(int);
static const int FLT_BIT = CHAR_BIT * sizeof(float);
static const int DBL_BIT = CHAR_BIT * sizeof(double);

//#ifndef bool
//typedef enum bool_t { false=0, true=1 } bool;
//#endif
typedef struct waveletFilter  * waveFilterPtr;
typedef struct wavePolyFilter * wavePolyPtr;
typedef struct waveLiftFilter * waveLiftPtr;
typedef struct dwtInfo *dwtInfPtr;

/* -------------- macro ----------------*/
#define NL printf("\n")
#define SWAP(A,B)    { typeof(A) tmp = A; A = B; B = tmp; }
#define XYN_MIN(x,y) ( (y<x) ? y : x )		// minimum
#define XYN_MAX(x,y) ( (y>x) ? y : x )		// maximum
#define XYN_DIV2(x)	 ( (x)/2 + (x)%2 )		// upper div by 2

/* Error Macros */
#define OK 		EXIT_SUCCESS
#define FAIL 	EXIT_FAILURE

#define XynFree(p) if ((p)!= NULL) {free(p);}

#define ErrPrintMsg(msg) {fprintf(stderr, "Error %s in %s Line %d, function %s\n",\
                            strerror(errno),__FILE__,__LINE__,__FUNCTION__);\
                  perror(msg);}
#define ErrPrint {fprintf(stderr, "Error %s in %s Line %d, function %s\n",\
                            strerror(errno),__FILE__,__LINE__,__FUNCTION__); }
#define ErrExitMsg(msg) { ErrPrintMsg(msg); errno=(errno==OK)?FAIL:errno; exit(errno); }
#define ErrExitMsgVA(format,...) { ErrPrint; fprintf(stderr, format, __VA_ARGS__); errno=(errno==OK)?FAIL:errno; exit (errno); }
#define ErrMsg(msg) { ErrPrintMsg(msg); errno=(errno==OK)?FAIL:errno; return(errno); }
#define ErrMsgVA(format,...) { ErrPrint; fprintf(stderr, format, __VA_ARGS__); errno=(errno==OK)?FAIL:errno; return(errno); }

/* profiling
#define showPartTime(x,t) showTotalTime("\t",x,t)
#define showTotalTime(tab,x,t) { printf("\n\t%s %-30s %lf\n",tab,x,t); fprintf(logFile,"\n\t%s %-30s %lf\n",tab,x,t); }
#define showTime(x) {double endTimeTmp=endTime(); showTotalTime("timing",x,endTimeTmp) }
*/
// define
//#define PI 3.14159265358979
//#define SQRT2 1.4142135623731

#define FILTER_LEN 20
#define POLY_LEN 10
#define PARSE_MAX_ARGUMENTS 1024
#define STRING_LEN 1024
#define DEF_BUFF_LEN 1024



typedef void* pointer;
typedef char xynString[STRING_LEN + 1];

/* xynString.c */
char *getLine ( FILE* fp );
char *skipWhiteSpaces ( FILE* fp );
char* skipWhiteLines ( FILE* fp );
char *getNewLine ( FILE* fp );
long getFileSize ( char* fname );

void stringMakeNull(xynString str);
int stringNull(const xynString str);
void convertToXynString(xynString strOut, const char *strIn);
void stringCopy(xynString str1, const xynString str2);
void stringSprintf(xynString str, const char *format, ...);

char* itoa(int value, char* result, int base);
int getIntFromFile ( FILE* fp, int n, int* out );
int getFloatFromFile ( FILE* fp, int n, float* out );


/*  xynParse.c  */
int parseParametersVA(int argc, char *argv[], const char *format, va_list ap);
int parseParameters(int argc, char *argv[], const char *format, ...);


#include "xynDisp.h"
#include "xynWave.h"
#include "xynRawVid.h"
#include "xynMath.h"

#undef _GNU_SOURCE
#endif /* LIBXYN_H_ */
