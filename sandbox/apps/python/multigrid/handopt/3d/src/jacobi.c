#include "smoothers.h"
#include <stdio.h>
#include <math.h>

//void jacobi( GRID &u, GRID &f, GRID &tmp, int nu ) {
long long int jacobi( GRID &u, GRID &f, GRID &tmp, int nu ) {
  int N = u.n;
  int lda = u.lda;
  int lda2 = u.lda * u.lda;
  double hh = u.h*u.h;
  double invhh = 1.0 / hh;
  double DinvXomega = hh/6.0 * 8.0/9.0;

  double* w = &(tmp.p[0][0][0]);
  double* a = &(u.p[0][0][0]);
  double* b = &(f.p[0][0][0]);
#ifdef USE_MM_ALLOC
  __assume_aligned(w,64);
  __assume_aligned(a,64);
  __assume_aligned(b,64);
#endif

long long int count = 0;

  for (int s=0; s<nu; s++) {
#pragma omp parallel for collapse(2)
    for (int i=1; i<=N; i++) {
      for (int j=1; j<=N; j++) {
	int lbv= i*lda2+j*lda+1;
	int ubv= i*lda2+j*lda+N;
	int ks = lbv-lda;
	int kn = lbv+lda;
	int ke = lbv-1;
	int kw = lbv+1;
	int kf = lbv-lda2;
	int kb = lbv+lda2;
	if (s%2==0) {
#pragma ivdep
#pragma vector always
	  for (int k=lbv; k<=ubv; k++) {
	    w[k] = a[k] - DinvXomega*((6.0*a[k]-a[kw]-a[ke]-a[kn]-a[ks]-a[kf]-a[kb])*invhh - b[k]);
	    kn++; ks++; ke++; kw++; kf++; kb++; 
        //count++;
	  }
	} else {
#pragma ivdep
#pragma vector always
	  for (int k=lbv; k<=ubv; k++) {
	    a[k] = w[k] - DinvXomega*((6.0*w[k]-w[kw]-w[ke]-w[kn]-w[ks]-w[kf]-w[kb])*invhh - b[k]);
	    kn++; ks++; ke++; kw++; kf++; kb++;
        //count++;
	  }
	}
      }
    }
  }
  if (nu%2==1) {
    double*** t = u.p;
    u.p = tmp.p;
    tmp.p = t;
  }

    return count;
}
