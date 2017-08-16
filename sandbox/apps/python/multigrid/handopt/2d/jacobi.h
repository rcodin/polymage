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

void jacobi( GRID* u, GRID* f, GRID* tmp, int nu ) {
  int N = u->n;
  int lda = u->lda;
  double hh = u->h*u->h;
  double invhh = 1.0 / hh;
  double DinvXomega = hh/4.0 * 4.0/5.0;

  //printf("invhh = %g\tDinvXomega = %g\n", invhh, DinvXomega);

  double* w = &(tmp->p[0][0]);
  double* a = &(u->p[0][0]);
  double* b = &(f->p[0][0]);
#ifdef USE_MM_ALLOC
  __assume_aligned(w,64);
  __assume_aligned(a,64);
  __assume_aligned(b,64);
#endif

  for (int k=0; k<nu; k++) {
#pragma omp parallel for// collapse(2)
    for (int i=1; i<=N; i++) {
      int lbv= i*lda+1;
      int ubv= i*lda+N;
      int ju = lbv+lda;
      int jb = lbv-lda;
      int jl = lbv-1;
      int jr = lbv+1;
      if (k%2==0) {
#pragma ivdep
#pragma vector always
	for (int j=lbv; j<=ubv; j++) {
	  w[j] = a[j] - DinvXomega*((4.0*a[j]-a[ju]-a[jl]-a[jb]-a[jr])*invhh - b[j]);
      //printf("%g %g %g %g %g %g %g\n", w[j], a[j], a[ju], a[jl], a[jb], a[jr], b[j]);
	  jr++; jl++; ju++; jb++;
	}
      } else {
#pragma ivdep
#pragma vector always
	for (int j=lbv; j<=ubv; j++) {
	  a[j] = w[j] - DinvXomega*((4.0*w[j]-w[ju]-w[jl]-w[jb]-w[jr])*invhh - b[j]);
	  jr++; jl++; ju++; jb++;
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
