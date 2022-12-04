/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 * 
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   Write data in (native-endian) SU format
 *------------------------------------------------------------------------*/

#include "write_su.h"
#include "su_gather.h"
#include "su_struct.h"
#include "logging.h"
#include <stdlib.h>


int su_write_trace(FILE *filep, const SUhead *header, const float *data)
{
  int ierr=0;
  size_t iwritten=0;
  unsigned short ns=0;
  
  if (!filep) {
    log_error("Cannot write SU data to disk, file pointer is NULL.\n");
    ierr = 1;
  }
  if ((!header) || (!data)) {
    log_error("Cannot write SU data to disk, buffers are NULL pointers.\n");
    ierr = 1;
  }
  ns = header->ns;
  if (ns<1) {
    log_error("Cannot write SU data to disk, trace length is zero.\n");
    ierr = 1;
  }

  if (ierr!=0) return ierr;

  iwritten = fwrite(header, sizeof(SUhead), 1, filep);
  if (iwritten != 1) {
    log_error("Short write while writing SU header to disk.\n");
    return 2;
  }

  iwritten = fwrite(data, sizeof(float), ns, filep);
  if (iwritten != ns) {
    log_error("Short write while writing SU data to disk.\n");
    return 3;
  }

  return ierr;
}

int su_write_gather(FILE *filep, const SUgather *gather)
{
  int ierr=0;
  size_t i;
  
  if (!gather) {
    log_error("Cannot write SU gather to disk, gather is NULL.\n");
    return 1;
  }
 
  for (i=0; i<gather->nt; ++i) {
    ierr += su_write_trace(filep, &(gather->header[i]), gather->data[i]);
    if (ierr!=0) break; // there is not much point in continuing
  }

  return ierr;
}




///////////////////////////////////////////////////////////////////////////////

#ifdef WRITE_SU_MAIN

#include <stdio.h>

// test main program //
int main() 
{
  FILE *fp = fopen("su_write_main_test.su", "wb");
  if (!fp) {
    printf("Could not open test SU file for writing.\n");
    exit(EXIT_FAILURE);
  }

  SUgather g;
  const size_t nt = 5;
  const unsigned short ns = 11;
  size_t count;

  malloc_SUgather(&g, nt, ns);
  
  for (size_t i=0; i<nt; ++i) {
    g.header[i].ns = ns;
    g.header[i].dt = 1000;
    g.header[i].ep = 1;
    g.header[i].tracl = i+1;
    g.header[i].tracf = i+1;
    g.header[i].ntr = (int)nt;
    g.header[i].sx = 42;
    g.header[i].gx = 40+i;
    g.header[i].offset = g.header[i].sx-g.header[i].gx;
    g.header[i].d1 = 0.001;
    count = 0;
    for (size_t j=0; j<ns; ++j) {
      g.data[i][j] = i*(ns+1)+count;
      ++count;
    }
  }

  int ierr = su_write_gather(fp, &g);
  if (ierr!=0) {
    printf("Error while writing gather to disk.\n");
    exit(EXIT_FAILURE);
  }

  free_SUgather(&g);
  
  if (fp) fclose(fp);

  exit(EXIT_SUCCESS);
}

#endif

///////////////////////////////////////////////////////////////////////////////




