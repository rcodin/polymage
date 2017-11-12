#ifndef intergrid_h_
#define intergrid_h_

#include "utils.h"
#include "grid.h"

/*
void restriction( GRID uf, GRID uc );
void prolong( GRID uf, GRID uc );
void correct( GRID u, GRID c );
*/
long long int restriction( GRID uf, GRID uc );
long long int prolong( GRID uf, GRID uc );
long long int correct( GRID u, GRID c );

#endif // intergrid_h_
