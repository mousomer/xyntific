/* string manipulation */

#include "libxyn.h"


char* skipWhiteSpaces ( FILE* fp )
{
	assert(fp);
    char* line=NULL;
    while ( 1 )
    {
    	if(line!=NULL)
    		free(line);
        line = getLine ( fp );
        if ( line )
        {
            if ( line[0]=='#' || line[0]=='\n')
            	continue;
            else
            	line = strtok ( line, " \t" );
        }
        else
            if ( ! feof ( fp ) )
            	continue;

        break;
    }
    return line;
}

char* skipWhiteLines ( FILE* fp )
{
	assert(fp);
    char* line=NULL;
    while ( 1 )
    {
    	if(line!=NULL)
    		free(line);

        line = getLine ( fp );
        if ( line )
        {
            if ( line[0]=='#' || line[0]=='\n')
            	continue;
        }
        else
            if ( ! feof ( fp ) )
            	continue;

        break;
    }
    return line;
}

char* getLine ( FILE* fp )
{
    assert(fp);
    static char* line = NULL;
    size_t len = 0;
    int nl = getline ( &line, &len, fp );
    assert ( nl != -1);
    assert( nl <= ( int ) len );

    if ( nl == EOF && feof ( fp ) )	return NULL;
    assert( nl != EOF );

    return line;
}

char* getNewLine ( FILE* fp )
{
	assert(fp);
	char* lineInp = getLine (fp);
	char* lineOut = (char*) malloc ( sizeof(char) * strlen(lineInp) );
	strcpy (lineOut, lineInp);
	return (lineOut);
}


int readIntText ( FILE* fp, long *n )
{
	assert(fp);
    char* line = getLine ( fp );
    sscanf ( line, "%ld", n );

    return (OK);
}

int readDblText ( FILE* fp, double *n )
{
	assert(fp);
    char* line = getLine ( fp );
    sscanf ( line, "%lf", n );

    return (OK);
}

int readNumberFromText ( FILE* fp, void* n, int nBytes, int len, int (*aToNum)(char*, void*) )
{
  assert( (fp!=NULL) && (n!=NULL) && (nBytes>0) && (len>=0) );
    char *line = getLine ( fp );
    char *num, *l = line;
    int okCode = 0;
    while ( len-- > 0 ) 
    {
        num = strtok ( l," \n\t" );
        assert( num );
        okCode = aToNum( num, n );
        assert (okCode);
        n += nBytes;
    }

    free (line);
    return OK;
}

int fileTrunc ( char* src, char* trg, long nBytes )
{
    assert( (trg!=NULL) && (src!=NULL) && (nBytes>=0) );

	FILE* fpInp = fopen ( src, "r" );
	FILE* fpOut = fopen ( trg , "w" );
    assert ( (fpOut!=NULL) && (fpInp!=NULL) );

    char c=0;

    while ( (nBytes-- > 0) && ((c=fgetc(fpInp)) != EOF) )
      putc( c, fpOut );

	fclose( fpInp );
	fclose( fpOut );
	return OK;
}

long getFileSize ( char* fName )
{
	struct stat buf;
    assert ( stat ( fName,&buf ) == 0 );
    return (buf.st_size);
}



void stringMakeNull(xynString str)
{
  if (str == NULL)
    return;

  str[0] = '\0';
}


int stringNull(const xynString str)
{
  if (str == NULL)
    return(1);
  return(strlen(str) <= 0);
}


void convertToXynString(xynString strOut, const char *strIn)
{
  if ((strIn == NULL) || (strOut == NULL))
    return;

  strncpy((char *)strOut, strIn, STRING_LEN);
  strOut[STRING_LEN] = '\0';
}


void stringCopy(xynString str1, const xynString str2)
{
  if (str1 == NULL)
    return;
  if (str2 == NULL)
    return;

  strncpy((char *)str1, (char *)str2, STRING_LEN);
  str1[STRING_LEN] = '\0';

}


void stringSprintf(xynString str, const char *format, ...)
{
  va_list ap;

  va_start(ap, format);


  vsprintf(str, format, ap);

//  vsnprintf(str, STRING_LEN, format, ap);

  str[STRING_LEN] = '\0';

  va_end(ap);
}


char* itoa(int value, char* result, int base)
{
		// check that the base if valid
		if (base < 2 || base > 36) { *result = '\0'; return result; }

		char* ptr = result, *ptr1 = result, tmp_char;
		int tmp_value;

		do {
			tmp_value = value;
			value /= base;
			*ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
		} while ( value );

		// Apply negative sign
		if (tmp_value < 0) *ptr++ = '-';
		*ptr-- = '\0';
		while(ptr1 < ptr) {
			tmp_char = *ptr;
			*ptr--= *ptr1;
			*ptr1++ = tmp_char;
		}
		return result;
	}



int getIntFromFile ( FILE* fp, int n, int* out )
{
	assert( fp!=NULL && out!=NULL );

	char* line = skipWhiteLines( fp );
	assert (line!=NULL);
	XYN_DPRINTF(DEBUG_HEADER,"%s line pointer %p, trying to fetch %d integers from line:\n\t%s\n", __FUNCTION__, line, n, line);

	char* lineCpy = line;
	while (n-->0)
	{
		*out = (int) strtol(lineCpy, &lineCpy, 10);
		XYN_DPRINTF(DEBUG_DETAIL,"got %d from %p. Further line: %s\n", *out, lineCpy, lineCpy );
		out++;
	}

	if (line != NULL)
		free (line);

	return OK;
}


int getFloatFromFile ( FILE* fp, int n, float* out )
{
	assert( fp!=NULL && out!=NULL );

	char* line = skipWhiteLines( fp );
	assert (line!=NULL);
	XYN_DPRINTF(DEBUG_HEADER,"%s line pointer %p (into %p), trying to fetch %d integers from line:\n\t%s\n",
			__FUNCTION__, line, out, n, line);

	char* lineCpy = line;
	while (n-->0)
	{
		*out = (float) strtof(lineCpy, &lineCpy);
		XYN_DPRINTF(DEBUG_DETAIL,"got %f from %p. Further line: %s\n", *out, lineCpy, lineCpy );
		out ++;
	}

	if (line != NULL)
		free (line);

	return OK;
}

