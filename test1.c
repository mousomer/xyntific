/*
 *
 * QccPack: Quantization, compression, and coding utilities
 * Copyright (C) 1997-2009  James E. Fowler
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */
#include "libxyn.h"

#define USG_STRING "[-prec %d:precision] %d:digits %s:type %f:fnum %d:dnum"

int precision=0, digits=0, dnum=0;
float fnum=0.0;
xynString type;

int main(int argc, char *argv[])
{
  if (parseParameters(argc, argv,
                         USG_STRING,
                         &precision,
                         &digits,
                         type,
                         &fnum,
                         &dnum))  return(FAIL);

	xynString fmt;

	if ( precision )
	{
		//stringSprintf(fmt, " /%%d.%d%c ", digits, precision, type[0] );
		printf("digits %d precision %d type %c\n format %s\n", digits, precision, type[0], fmt );
		printf(" %*.*f ", digits, precision, fnum);
	}
	else
	{
		//stringSprintf(fmt, " /%%d%c ", digits, type[0] );
		printf("digits %d type %c\n format %s\n", digits, type[0], fmt );
		printf(" %*d ", digits, dnum);
	}

	printf("%s\n", fmt);

  
  

  return OK;
}
