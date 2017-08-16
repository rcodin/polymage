#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define pmax(x,y)    ((x) > (y)? (x) : (y))
#define pmin(x,y)    ((x) < (y)? (x) : (y))

#include "smoothers.h"
#include <math.h>

// cannot be 8 or below
#ifndef TS
#define TS 16
#endif

#ifndef T3
#define T3 64
#endif

#ifndef T4
#define T4 256
#endif


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

    long long int count = 0;

#ifdef USE_MM_ALLOC
  __assume_aligned(w,64);
  __assume_aligned(a,64);
  __assume_aligned(b,64);
#endif

  if ((N >= 1) && (nu >= 1)) {
    for (int t1=-1;t1<=floord(nu-1,TS);t1++) {
      int lbp=pmax(ceild(t1,2),ceild(TS*t1-nu+2,TS));
      int ubp=pmin(floord(nu+N-1,TS),floord((TS/2)*t1+N+(TS/2)-1,TS));
#pragma omp parallel for
      for (int t2=lbp;t2<=ubp;t2++) {
	for (int t3=pmax(pmax(0,ceild(t1-1,2)),ceild(TS*t2-N-(T3-2),T3));t3<=pmin(pmin(floord(nu+N-1,T3),floord((TS/2)*t1+N+T3-1,T3)),floord(TS*t2+N+T3-2,T3));t3++) {
	  for (int t4=pmax(pmax(pmax(0,ceild(t1-((T4/4)-1),(T4/4))),ceild(TS*t2-N-(T4-2),T4)),ceild(T3*t3-N-(T4-2),T4));t4<=pmin(pmin(pmin(floord(nu+N-1,T4),floord((TS/2)*t1+N+TS-1,T4)),floord(TS*t2+N+TS-2,T4)),floord(T3*t3+N+T3-2,T4));t4++) {

	    for (int t5=pmax(pmax(pmax(pmax(pmax(0,(TS/2)*t1),TS*t2-N),T3*t3-N),T4*t4-N),TS*t1-TS*t2+1);t5<=pmin(pmin(pmin(pmin(pmin(nu-1,(TS/2)*t1+TS-1),TS*t2+TS-2),T3*t3+T3-2),T4*t4+T4-2),TS*t1-TS*t2+N+(TS-1));t5++) {
#pragma loop_count min(1),max(8),avg(4)
	      for (int t6=pmax(pmax(TS*t2,t5+1),-TS*t1+TS*t2+2*t5-(TS-1));t6<=pmin(pmin(TS*t2+TS-1,t5+N),-TS*t1+TS*t2+2*t5);t6++) {
#pragma loop_count min(1),max(8),avg(4)
		for (int t7=pmax(T3*t3,t5+1);t7<=pmin(T3*t3+T3-1,t5+N);t7++) {
		  int lbv= (-t5+t6)*lda2 + (-t5+t7)*lda -t5+pmax(T4*t4,t5+1);
		  int ubv= (-t5+t6)*lda2 + (-t5+t7)*lda -t5+pmin(T4*t4+T4-1,t5+N);
		  int ks = lbv-lda;
		  int kn = lbv+lda;
		  int ke = lbv-1;
		  int kw = lbv+1;
		  int kf = lbv-lda2;
		  int kb = lbv+lda2;
		  if (t5%2==0) {
#pragma loop_count min(1),max(64),avg(32)
#pragma ivdep
#pragma vector always
		    for (int k=lbv;k<=ubv;k++) {
		      w[k] = a[k] - DinvXomega*((6.0*a[k]-a[kw]-a[ke]-a[kn]-a[ks]-a[kf]-a[kb])*invhh - b[k]);
		      kn++; ks++; ke++; kw++; kf++; kb++; 

    //count++;
		    }
		  } else {
#pragma loop_count min(1),max(64),avg(32)
#pragma ivdep
#pragma vector always
		    for (int k=lbv;k<=ubv;k++) {
		      a[k] = w[k] - DinvXomega*((6.0*w[k]-w[kw]-w[ke]-w[kn]-w[ks]-w[kf]-w[kb])*invhh - b[k]);
		      kn++; ks++; ke++; kw++; kf++; kb++;

   // count++;
		    }
		  }
		}
	      }
	    }
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
