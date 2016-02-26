/*
 * xynString.h
 *
 *  Created on: Jan 27, 2011
 *      Author: mousomer
 */

#ifndef XYNSTRING_H_
#define XYNSTRING_H_


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define DEF_BUFF_LEN 1024

char *getLine ( FILE* fp );
char *skipWhiteSpaces ( FILE* fp );
char *getNewLine ( FILE* fp );
long getFileSize ( char* fname );

int getIntFromFile ( FILE* fp, int n, int* out );

#endif /* XYNSTRING_H_ */
