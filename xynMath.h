/*
 * xynmath.h
 *
 *  Created on: 07/01/2011
 *      Author: lilah215
 */

#ifndef XYNMATH_H_
#define XYNMATH_H_

unsigned long pow2_i( unsigned int m );
int 	log2_l		( unsigned long m );
int 	log2whole_l	( unsigned long m );
long 	powNM		( int n, int m );
double 	pow_dbyi 	( double f, int n );

/** integer vector operations **/
long 	sum_i 				( int* v, long n );
long 	prod_i 				( int* x, long n );
int 	addVectScal_i		( int* v, int scalar, long n );
int 	addVectScalWin_i	( int* v, int scalar, int nC, int winC, int winR, int hOfs, int vOfs );
int 	mulVectScal_i		( int* v, int scalar, long n );
int 	mulVectScalJump_i	( int* v, int scalar, long n, int jmp );
int 	mulVectScal_f		( float* v, float scalar, long n );
int 	mulVectScalJump_f	( float* v, float scalar, long n, int jmp );
int 	mulVectScalJump_d	( double* v, double scalar, long n, int jmp );

long	dotProd_i		( int* v, int* u, long len );
long 	dotProdJump_i 	( int* v, int jv, int* u, int ju, long len );
int 	cpyJump_i 		( int* inp, int ji, int* out, int jo, long len );


/** integer vector operations **/
int		addVectorsTriple_d	( double* x0, double* xp, double* xn, int len, double cp, double cn );
int		addVectorsTriple_f	( float* x0, float* xp, float* xn, int len, float cp, float cn );
int		addVectorsTriple_i	( int* x0, int* xp, int* xn, int len, int cp, int cn, int downScaler );
int		addFramesTriple_i ( int* fr0, int n0, int* frP, int np, int* frN, int nN, segmInfPtr sp,
								register int cp, register int cn, int downScaler );
double	dotProd_f 		( float* v, float* u, long len );
double	dotProdJump_f 	( float* v, int jv, float* u, int ju, long len );
int 	cpyJump_f 		( float* inp, int ji, float* out, int jo, long len );

/** lower level comp **/
int 	lowLen	 	( long length, int level );
int* 	SBlengths 	( long length, int level );

/* moments and central moment. Energy is 2nd (non central) moment. */
double 	mean_i 		( int* y, long n );
double 	mean_f 		( float* y, long n );
double 	meanDiv_i 	( int* y, long n );
double 	meanDiv_f 	( float* y, long n );
double 	moment_i	( int* y, long n, int moment, int central );
double 	moment_f	( float* y, long n, int moment, int central );
double 	momentDiv_i	( int* y, long n, int moment, int central );
double 	momentDiv_f	( float* y, long n, int moment, int central );
double 	meanWin_i 	( int* y, int nC, int winC, int winR, int hOfs, int vOfs );
double 	momentWin_i ( int* y, int moment, int central, int nC, int winC, int winR, int hOfs, int vOfs );



double 	diff_i 		( int* x, int* y, long n, double(*diffunc)(int,int) );
double 	diffWin_i 	( int* x, int* y, int nC, int winC, int winR, int hOfsX, int vOfsX,
					  int hOfsY, int vOfsY, double(*diffunc)(int,int) );
/*
int filterQMF_i		(int* inp, int* out, int len, int zeroPos);
int filterQMF_f		(float* inp, float* out, int len, int zeroPos);
*/

int linearVect_f 	(float* data, long len, float min, float max );
int constVect_f 	(float* data, long len, float min, float max);
int powerVect_f		(float* data, long len, float min, float max, int n );

int 	randVector_i	( int* ip, int n );
int*    genRandVector_i	( long len, int max);
float* 	genRandVector_f ( long len, float max );

/* conversions */
int convert_ftoi ( int* intData, const float* floatData, long len );
int convert_ftog ( double* dblData, const float* floatData, long len );


#endif /* XYNMATH_H_ */
