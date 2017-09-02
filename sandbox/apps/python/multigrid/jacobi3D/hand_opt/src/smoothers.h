#ifndef smoothers_h_
#define smoothers_h_

#include "utils.h"
#include "grid.h"


//void jacobi( GRID &u, GRID &f, GRID &tmp, int nu);
long long int jacobi( GRID &u, GRID &f, GRID &tmp, int nu);
//void chebyshev( GRID &u, GRID &f, GRID &tmp, int nu, double lmin, double lmax);

#endif // smoothers_h_
