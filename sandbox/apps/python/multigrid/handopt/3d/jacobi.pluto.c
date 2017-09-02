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
  for (int t1=-1;t1<=floord(nu-2,8);t1++) {
    int lbp=pmax(ceild(t1,2),ceild(16*t1-nu+3,16));
    int ubp=pmin(floord(nu+N-2,16),floord(8*t1+N+7,16));
#pragma omp parallel for
    for (int t2=lbp;t2<=ubp;t2++) {
      for (int t3=pmax(pmax(0,ceild(t1-1,2)),ceild(16*t2-N-14,16));t3<=pmin(pmin(floord(nu+N-2,16),floord(8*t1+N+15,16)),floord(16*t2+N+14,16));t3++) {
        for (int t4=pmax(pmax(pmax(0,ceild(t1-124,125)),ceild(16*t2-N-998,1000)),ceild(16*t3-N-998,1000));t4<=pmin(pmin(pmin(pmin(floord(8*t1-8*t2+N+7,500),floord(nu+N-2,1000)),floord(8*t1+N+15,1000)),floord(16*t2+N+14,1000)),floord(16*t3+N+14,1000));t4++) {
          for (int t5=pmax(pmax(pmax(pmax(pmax(0,8*t1),16*t2-N),16*t3-N),1000*t4-N),16*t1-16*t2+1);t5<=pmin(pmin(pmin(pmin(pmin(nu-2,8*t1+15),16*t2+14),16*t3+14),1000*t4+998),16*t1-16*t2+N+15);t5++) {
            for (int t6=pmax(pmax(16*t2,t5+1),-16*t1+16*t2+2*t5-15);t6<=pmin(pmin(16*t2+15,t5+N),-16*t1+16*t2+2*t5);t6++) {
              for (int t7=pmax(16*t3,t5+1);t7<=pmin(16*t3+15,t5+N);t7++) {
		  int lbv= (-t5+t6)*lda2 + (-t5+t7)*lda -t5+pmax(1000*t4,t5+1);
		  int ubv= (-t5+t6)*lda2 + (-t5+t7)*lda -t5+pmin(1000*t4+999,t5+N);
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
/*
                lbv=pmax(1000*t4,t5+1);
                ubv=pmin(1000*t4+999,t5+N);
#pragma ivdep
#pragma vector always
                for (int t8=lbv;t8<=ubv;t8++) {
                  a[( t5 + 1) % 2][ (-t5+t6)][ (-t5+t7)][ (-t5+t8)] = (a[ t5 % 2][ (-t5+t6)][ (-t5+t7)][ (-t5+t8)] - (c * (((((((((6.0 * a[ t5][ (-t5+t6)][ (-t5+t7)][ (-t5+t8)]) - a[ t5 % 2][ (-t5+t6) - 1][ (-t5+t7)][ (-t5+t8)]) - a[ t5 % 2][ (-t5+t6) + 1][ (-t5+t7)][ (-t5+t8)]) - a[ t5 % 2][ (-t5+t6)][ (-t5+t7) - 1][ (-t5+t8)]) - a[ t5 % 2][ (-t5+t6)][ (-t5+t7) + 1][ (-t5+t8)]) - a[ t5 % 2][ (-t5+t6)][ (-t5+t7)][ (-t5+t8) - 1]) - a[ t5 % 2][ (-t5+t6)][ (-t5+t7)][ (-t5+t8) + 1]) * invhh) - b[ t5 % 2][ (-t5+t6)][ (-t5+t7)][ (-t5+t8)])));;
                }*/
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
