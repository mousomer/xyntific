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

#define USG_STRING "%s:fName"
xynString fName;

int main(int argc, char *argv[])
{
  if (parseParameters(argc, argv, USG_STRING,
                         &fName))  return(FAIL);

  segmInfPtr sp = genSegmentInfo (fName, 1);
  assert(sp);
  XYN_DPRINTF(DEBUG_PROCESS, "generated segment struct pointer %p\n", sp );
  printSegmentInfo ( stdout, sp );
  
  free(sp);

  return OK;
}
