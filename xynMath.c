/**
*  @file xynMath.c
*
*		The file contains basic mathematical functions for integer data.
*
* @author: Omer Moussaffi <mousomer@gmail.com>, (C) 2011
*
* Copyright: See COPYING file that comes with this distribution
*
*/
#include "libxyn.h"

//static int addFramesTriple_i ( int* fr0, int n0, int* frP, int np, int* frN, int nN, segmInfPtr sp, register int cp, register int cn, int downScaler );


/* long integer powers and logarithms */
unsigned long pow2_i ( unsigned int m )
{
	assert ( m < LNG_BIT );
    return 1 << m;
}

int log2_l ( unsigned long m )
{
	int k = 0;
	while ( m >= 1 )
	{
		k++;
		m = m >> 1;
	}
    return k;
}

int log2whole_l ( unsigned long m )
{
	int k = 0;
	long n = 1;
	while ( n < m )
	{
		n = n << 1;
		k++;
	}
	return k;
}

/* return n^m */
long powNM ( int n, int m )
{
    long k = 1;
    while ( m-- > 0 )
        k *= n;
    return k;
}

/* return f^n for double and integer */
double pow_dbyi ( double f, int n )
{
	assert(n>=0);
	n--;
	double out = f;
	while (n-- > 0)
		out *= f;
	return out;
}


/** vector operations **/
/* sum vector elements */
long sum_i ( int* v, long n )
{
    assert( v && n>0 );

	long sum = 0;
	for ( long i=0; i<n; i++ )
		sum += *v++;

	return sum;
}

/* multiply vector elements */
long prod_i ( int* x, long n )
{
    assert( x && n>0 );
    long m = *x;
    for ( int i=1; i<n; i++ )
        m *= *x++;

    return m;
}

/* add scalar to vector elements */
int addVectScal_i ( int* v, int scalar, long n )
{
    assert( v && n>0 );

    	for ( long i=0; i<n; i++ )
        *v++ += scalar;

    return OK;
}

/* add scalar to a window of vector elements */
int addVectScalWin_i ( int* v, int scalar, int nC, int winC, int winR, int hOfs, int vOfs )
{
    assert( v && nC>0 && winC>=0 && winR>=0 && hOfs>70 && vOfs>=0 );

	v += nC*vOfs+hOfs;
	int res = nC - winC;

	for ( int j=0; j<winR; j++, v+=res )
		for ( int i=0; i<winC; i++, v++ )
			*v += scalar;

	return OK;
}

/* multiply vector elements by scalar */
int mulVectScal_i ( int* v, int scalar, long n )
{
    assert( v && n>0 );

    for ( long i=0; i<n; i++ )
        *v++ *= scalar;

    return OK;
}

/* multiply vector elements by scalar with a jump*/
int mulVectScalJump_i ( int* v, int scalar, long n, int jmp )
{
    assert( v && n>0 );

    for ( long i=0; i<n; i+=jmp )
    {
        *v *= scalar;
        v += jmp;
    }

    return OK;
}

/* multiply vector elements by scalar */
int mulVectScal_f ( float* v, float scalar, long n )
{
    assert( v && n>0 );

    for ( long i=0; i<n; i++ )
        *v++ *= scalar;

    return OK;
}

/* multiply vector elements by scalar with a jump*/
int mulVectScalJump_f ( float* v, float scalar, long n, int jmp )
{
    assert( v && n>0 );

    for ( long i=0; i<n; i+=jmp )
    {
        *v *= scalar;
        v += jmp;
    }

    return OK;
}

int mulVectScalJump_d ( double* v, double scalar, long n, int jmp )
{
    assert( v && n>0 );

    for ( long i=0; i<n; i+=jmp )
    {
        *v *= scalar;
        v += jmp;
    }

    return OK;
}

/** dot products: _i and _f **/
long dotProd_i 	( int* v, int* u, long len )
{
	assert (v && u && len>0);
	long product = 0;

	while (len-- > 0)
		product += *v++ * *u++;

	return product;
}

long dotProdJump_i 	( int* v, int jv, int* u, int ju, long len )
{
	assert (v && u && len>0);
	long product = 0;

	while (len-- > 0)
	{
		product += *v * *u;
		v+=jv;
		u+=ju;
	}

	return product;
}

double dotProd_f 	( float* v, float* u, long len )
{
	assert (v && u && len>0);
	double product = 0;

	while (len-- > 0)
		product += *v++ * *u++;

	return product;
}

double dotProdJump_f 	( float* v, int jv, float* u, int ju, long len )
{
	assert (v && u && len>0);
	double product = 0;

	while (len-- > 0)
	{
		product += *v * *u;
		v+=jv;
		u+=ju;
	}

	return product;
}

int cpyJump_f 	( float* inp, int ji, float* out, int jo, long len )
{
	assert (inp!=NULL && out!=NULL && len>0);

	while (len-- > 0)
	{
		*out = *inp;
		inp += ji;
		out += jo;
	}

	return OK;
}

int cpyJump_i 	( int* inp, int ji, int* out, int jo, long len )
{
	assert (inp!=NULL && out!=NULL && len>0);

	while (len-- > 0)
	{
		*out = *inp;
		inp += ji;
		out += jo;
	}

	return OK;
}

/** vector operations **/
int		addVectorsTriple_f	( float* x0, float* xp, float* xn, int len, float cp, float cn )
{
	assert (x0!=NULL && xp!=NULL && xn!=NULL );
	while (len-->0)
		*x0++ += cp*(*xp++) + cn*(*xn++);

	return OK;
}

int		addVectorsTriple_d	( double* x0, double* xp, double* xn, int len, double cp, double cn )
{
	assert (x0!=NULL && xp!=NULL && xn!=NULL );
	while (len-->0)
		*x0++ += cp*(*xp++) + cn*(*xn++);

	return OK;
}


int		addVectorsTriple_i	( int* x0, int* xp, int* xn, int len, int cp, int cn, int downScaler )
{
	assert (x0!=NULL && xp!=NULL && xn!=NULL );
	while (len-->0)
		*x0++ +=  ( cp*(*xp++) + cn*(*xn++) ) / downScaler;

	return OK;
}

int		addFramesTriple_i ( int* fr0, int n0, int* frP, int np, int* frN, int nN, segmInfPtr sp, register int cp, register int cn, int downScaler )
{
	assert (fr0!=NULL && frP!=NULL && frN!=NULL && sp!=NULL );
	//long len = (sp->frameSize)[sp->yuvCurrent];


	long len = sp->macroBlockYUV;
	while (len-->0)
	{
		*fr0++ +=  ( cp*(*frP++) + cn*(*frN++) ) / downScaler;
	}

	return OK;
}

/** subband length calculators **/
int lowLen ( long length, int level )
{
	for ( int i=0; i<level; i++)
		length = (length+1)/2;

	return length;
}

int* SBlengths ( long length, int level )
{
    assert( level>=0 && length>0 );

	int* lens = (int*) malloc ( sizeof(int) * (level+1) );
	assert( lens );

	for ( int i=level; i>0; i-- )
	{
		lens[i] = length/2;
		length = (length+1)/2;
	}
	lens[0] = length;

	return lens;
}

/* moments and central moment. Energy is 2nd (non central) moment.
 * versions: 	int_i
 * 				float_f */
double mean_i ( int* y, long n )
{
    assert( y && n>0 );

	double sum = 0.0;
	for ( long i=0; i<n; i++ )
		sum += *y++;
	sum /= n;

	return sum;
}

double mean_f ( float* y, long n )
{
    assert( y && n>0 );

	double sum = 0.0;
	for ( long i=0; i<n; i++ )
		sum += *y++;
	sum /= n;

	return sum;
}

double meanDiv_i ( int* y, long n )
{
    assert( y && n>0 );

	double sum = 0.0;
	for ( long i=0; i<n; i++ )
		sum += (double) *y++ / n;

	return sum;
}

double meanDiv_f ( float* y, long n )
{
    assert( y && n>0 );

	double sum = 0.0;
	for ( long i=0; i<n; i++ )
		sum += (double) *y++ / n;

	return sum;
}

double moment_i( int* y, long n, int moment, int central )
{
    assert( y && n>0 );

	double mean = 0.0, sum = 0.0;
	if (central)
	   	mean = mean_i ( y, n );

	for ( long i=0; i<n; i++,y++ )
		sum += (double) powNM( *y - mean, moment );
	sum /= n;

	return sum;
}

double moment_f( float* y, long n, int moment, int central )
{
    assert( y && n>0 );

	double mean = 0.0, sum = 0.0;
	if (central)
	   	mean = mean_f ( y, n );

	for ( long i=0; i<n; i++,y++ )
		sum += (double) pow( *y - mean, moment );
	sum /= n;

	return sum;
}

double momentDiv_i( int* y, long n, int moment, int central )
{
	assert( y && n>0 );

	double sum = 0.0, mean = 0.0;
	if (central)
	   	mean = meanDiv_i ( y, n );

	for ( long i=0; i<n; i++,y++ )
		sum += (double) powNM( *y - mean, moment )/n;

	return sum;
}

double momentDiv_f( float* y, long n, int moment, int central )
{
	assert( y && n>0 );

	double sum = 0.0, mean = 0.0;
	if (central)
	   	mean = meanDiv_f ( y, n );

	for ( long i=0; i<n; i++,y++ )
		sum += (double) pow( *y - mean, moment )/n;

	return sum;
}

double meanWin_i ( int* y, int nC, int winC, int winR, int hOfs, int vOfs )
{
	assert( y && nC>0 && winC>=0 && winR>=0 );

	double sum = 0.0;
	y += nC*vOfs+hOfs;
	int res = nC - winC;

	for ( int j=0; j<winR; j++, y+=res )
		for ( int i=0; i<winC; i++, y++ )
			sum += *y;
	sum /= winC*winR;

	return sum;
}

double momentWin_i ( int* y, int moment, int central, int nC, int winC, int winR, int hOfs, int vOfs )
{
	assert( y && nC>0 && winC>=0 && winR>=0 );

    double mean = 0.0;
    if (central)
    	mean = meanWin_i ( y, nC, winC, winR, hOfs, vOfs );

	double sum = 0.0;
	y+=nC*vOfs+hOfs;
	int res = nC - winC;

	for ( int j=0; j<winR; j++, y+=res )
		for ( int i=0; i<winC; i++, y++ )
			sum += powNM( *y - mean, moment );
	sum /= winC*winR;

	return sum;
}



/** differentiators **/
double diff_i ( int* x, int* y, long n, double(*diffunc)(int,int) )
{
	assert( x && y && n>0 );

	double sum = 0;
	for ( long i=0; i<n; i++,x++,y++ )
		sum += (*diffunc)( *x, *y );
	sum /= n;

	return sum;
}

double diffWin_i ( int* x, int* y, int nC, int winC, int winR, int hOfsX, int vOfsX,
				   int hOfsY, int vOfsY, double(*diffunc)(int,int) )
{
	assert( x && y && nC>0 && winC>=0 && winR>=0 );
	double sum = 0;

	x+=nC*vOfsX+hOfsX;
	y+=nC*vOfsY+hOfsY;

	int res = nC - winC;

	for ( int j=0; j<winR; j++, x+=res, y+=res )
		for ( int i=0; i<winC; i++, x++, y++ )
			sum += (*diffunc)( *x, *y );

	sum /= winC*winR;

	return sum;
}

int diffAbs(int a, int b)
{
	return abs(a-b);
}
int diffSqr(int a, int b)
{
	static int c;
	c = a-b;
	return (c*c);
}
int diffSub(int a, int b)
{
	return (a-b);
}
int diffEnt(int a, int b)
{
	static int nBits, diff;
	diff = abs(a-b); nBits = 0;
	while(diff)
	{
		nBits += (diff % 2);
		diff = diff>>1;
	}
	return nBits;
}

/** vector generators **/

/* linear vector */
int linearVect_f (float* data, long len, float min, float max)
{
	*data++ = min;
	len--;

	double increment = ( (double) max - min ) / (double) len;
	double a = (double) min + increment;

	while ( len-->0 )
	{
		*data++ = (float) a;
		a += increment;
	}
	return OK;
}

/* constant vector */
int constVect_f(float* data, long len, float min, float max)
{
	min = (max-min)/2;

	while ( len-->0 )
		*data++ = min;
	return OK;
}

/**
 * y(0) = min
 * y(len-1) = max = min + a*(len-1)^n
 * -> p = (max-min)/(len-1)^n   */
int powerVect_f(float* data, long len, float min, float max, int n )
{
	*data++ = min; len--;
	double increment = (double) (max - min ) * pow(len,-n);
	XYN_DPRINTF(DEBUG_DETAIL, "%f -> %f\t. %f / %f = %e\n", min, max, max-min, pow(len,(double)n) );
	XYN_DPRINTF(DEBUG_HEADER, "%s length %ld. xi=(%e)+(%e)*i^(%d)\n", __FUNCTION__, len+1, min, increment, n);

	double a = (double) min;

	for (int i = 1; i< len; i++ )
	{
		*data++ = (float) ( a + increment*pow_dbyi(i,n) );
		XYN_DPRINTF(DEBUG_DETAIL, " %f ", data[-1]);
	}

	*data = max;
	XYN_DPRINTF(DEBUG_DETAIL, " %f \n", data[0]);
	return OK;
}


/* random vector generators */
int randVector_i( int* ip, int n )
{
	assert( n>0 ); assert( ip );

	 srand((unsigned)(time(0)));
	 while ( n-- > 0 )
		 *ip++ = rand();

	 return (OK);
}

int* genRandVector_i ( long len, int max )
{
	srand( (unsigned) time(0) );
	int* data = (int*) malloc ( sizeof(int) * len );
	assert (data!=NULL);

	printf("RAND_MAX %d\n", RAND_MAX);
	int* cpy = data;
	while ( len-- > 0 )
	{
		*cpy++ = (int) ( (1.0/RAND_MAX) * rand () * max );
	}

	return data;
}

float* genRandVector_f ( long len, float max )
{
	srand( (unsigned) time(0) );
	float* data = (float*) malloc ( sizeof(float) * len );
	assert (data!=NULL);

	printf("RAND_MAX %d\n", RAND_MAX);
	float* cpy = data;
	while ( len-- > 0 )
	{
		*cpy++ = (float) ( (1.0/RAND_MAX) * rand () * max );
	}

	return data;
}

/* conversions */
int convert_ftoi ( int* intData, const float* floatData, long len )
{
	while (len-->0)
		*intData++ = lrintf (*floatData++);

	return OK;
}

int convert_ftog ( double* intData, const float* floatData, long len )
{
	while (len-->0)
		*intData++ = (double) (*floatData++);

	return OK;
}
