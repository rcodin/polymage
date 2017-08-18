/*
-----------------------------------------------------------------------
Copyright 2013 Pieter Ghysels, University of Antwerp

Contact: ghyselsp@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define pmax(x,y)    ((x) > (y)? (x) : (y))
#define pmin(x,y)    ((x) < (y)? (x) : (y))

// cannot be 8 or below
#ifndef TS
#define TS 32
#endif

#ifndef T3
#define T3 64
#endif

void jacobi( GRID* u, GRID* f, GRID* tmp, int nu ) {

  //printf("Tile sizes: %d %d %d\n", TS, TS, T3);

  int i, j, k;
  int N = u->n;
  int lda = u->lda;
  double hh = u->h*u->h;
  double invhh = 1.0 / hh;
  double DinvXomega = hh/4.0 * 4.0/5.0;

  double* w = &(tmp->p[0][0]);
  double* a = &(u->p[0][0]);
  double* b = &(f->p[0][0]);
#ifdef USE_MM_ALLOC
  __assume_aligned(w,64);
  __assume_aligned(a,64);
  __assume_aligned(b,64);
#endif

  if ((N >= 1) && (nu >= 1)) {
    for (int t1=-1;t1<=floord(nu-1,TS/2);t1++) {
      int lbp=pmax(ceild(t1,2),ceild(TS*t1-nu+2,TS));
      int ubp=pmin(floord(nu+N-1,TS),floord((TS/2)*t1+N+(TS/2-1),TS));
#pragma omp parallel for
      for (int t2=lbp;t2<=ubp;t2++) {
	for (int t3=pmax(pmax(0,ceild(t1-1,2)),ceild(TS*t2-N-(T3-2),T3));t3<=pmin(pmin(floord(nu+N-1,T3),floord((TS/2)*t1+N+(TS-1),T3)),floord(TS*t2+N+(TS-2),T3));t3++) {
	  for (int t4=pmax(pmax(pmax(pmax(0,(TS/2)*t1),TS*t2-N),T3*t3-N),TS*t1-TS*t2+1);t4<=pmin(pmin(pmin(pmin(nu-1,(TS/2)*t1+(TS-1)),TS*t2+(TS-2)),T3*t3+(T3-2)),TS*t1-TS*t2+N+(TS-1));t4++) {
#pragma loop_count min(1),max(TS),avg(TS/2)
	    for (int t5=pmax(pmax(TS*t2,t4+1),-TS*t1+TS*t2+2*t4-(TS-1));t5<=pmin(pmin(TS*t2+(TS-1),t4+N),-TS*t1+TS*t2+2*t4);t5++) {
	      int lbv=(-t4+t5)*lda-t4+pmax(T3*t3,t4+1);
	      int ubv=(-t4+t5)*lda-t4+pmin(T3*t3+(T3-1),t4+N);
	      int t6l=lbv-1;
	      int t6r=lbv+1;
	      int t6u=lbv+lda;
	      int t6b=lbv-lda;
	      if (t4%2==0) {
#pragma loop_count min(1),max(TS),avg((TS/2))
#pragma ivdep
#pragma vector always
		for (int t6=lbv;t6<=ubv;t6++) {
		  w[t6]=a[t6]-DinvXomega*((4.0*a[t6]-a[t6b]-a[t6l]-a[t6u]-a[t6r])*invhh-b[t6]);
		  t6l++; t6r++; t6u++; t6b++;
		}
	      } else {
#pragma loop_count min(1),max(TS),avg((TS/2))
#pragma ivdep
#pragma vector always
		for (int t6=lbv;t6<=ubv;t6++) {
		  a[t6]=w[t6]-DinvXomega*((4.0*w[t6]-w[t6b]-w[t6l]-w[t6u]-w[t6r])*invhh-b[t6]);
		  t6l++; t6r++; t6u++; t6b++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  if (nu%2==1) {
    double** t = u->p;
    u->p = tmp->p;
    tmp->p = t;
  }
}
