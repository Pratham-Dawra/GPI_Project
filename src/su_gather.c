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
 *   Handle collections of SU traces as part of a SUgather
 *------------------------------------------------------------------------*/

#include "su_gather.h"
#include "logging.h"
#include "util.h"
#include <stdlib.h>


void malloc_SUgather(SUgather *gather, size_t nt, unsigned short ns)
{
  // allocate header[nt]
  gather->header = (SUhead*)calloc(nt, sizeof(SUhead));
  if (!gather->header) log_fatal("Could not allocate trace header memory buffer for SUgather.\n");

  gather->data = (float**)malloc2d(nt, (size_t)ns, sizeof(float));
  gather->nt = nt;
  gather->ns = ns;
  
  return;
}

void free_SUgather(SUgather *gather)
{
  if ((0==gather->nt) || (0==gather->ns) || (!gather->header) || (!gather->data)) {
    log_warn("Attempt to free SUgather not previously allocated. Request ignored.\n");
    return;
  }
  
  free(gather->data); 
  free(gather->header);
  gather->data = NULL;
  gather->header = NULL;
  gather->nt = 0;
  gather->ns = 0;
  
  return;
}




///////////////////////////////////////////////////////////////////////////////

#ifdef SU_GATHER_MAIN

int main() 
{
  SUgather g;
  const size_t nt = 4;
  const unsigned short ns = 11;

  malloc_SUgather(&g, nt, ns);
  
  size_t count = 0;
  for (size_t i=0; i<nt; ++i) {
    for (size_t j=0; j<ns; ++j) {
      g.data[i][j] = count++;
    }
  }
  for (size_t i=0; i<nt; ++i) {
    for (size_t j=0; j<ns; ++j) {
      printf("trace: %2ld, sample: %2ld, address: %p, value: %f\n", i, j, &(g.data[i][j]), g.data[i][j]);
    }
    printf("\n");
  }

  free_SUgather(&g);

  exit(EXIT_SUCCESS);
}

#endif

///////////////////////////////////////////////////////////////////////////////
