/*
 * xynMem.h
 *
 *  Created on: Mar 11, 2011
 *      Author: mousomer
 */

#ifndef XYNMEM_H_
#define XYNMEM_H_

int xynConvert_c2i ( char* inp, float* out, long len );
int xynConvert_c2f ( char* inp, float* out, long len );
int xynConvert_i2c ( char* inp, float* out, long len );
int xynConvert_i2f ( char* inp, float* out, long len );
int xynConvert_f2c ( char* inp, float* out, long len );
int xynConvert_f2i ( char* inp, float* out, long len );

#endif /* XYNMEM_H_ */
