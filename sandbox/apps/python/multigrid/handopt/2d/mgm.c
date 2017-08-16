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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>


typedef enum {V,W,S} CYCLE;
typedef enum {JACOBI,CHEBYSHEV} METHOD;

typedef struct{ double **p; int n; double h; int lda; } GRID;
typedef GRID *GRIDLIST;

void jacobi( GRID* u, GRID* f, GRID* tmp, int nu );
#ifdef PLUTO
#include "jacobi.pluto.h"
#else
#include "jacobi.h"
#endif

double get_time();
void init_rhs( GRID u);
void init_border( GRID u );
void init_solution( GRID u );
void create_grid( GRID *grid, int n );
void free_grid( GRID *grid );
void create_grid_list( GRIDLIST *list, int n, int lmax );
void free_grid_list( GRIDLIST *list, int lmax );
void print_grid( GRID u, char *c );
void zero_grid( GRID u );
void init_grid( GRID u, double a );
void random_grid( GRID u );
double error( GRID u1, GRID u2 );
void correct( GRID u, GRID c );
void get_options( int *iters, int *nc, int *lmax, int *nu1, int *nu2, int *nuc, int *print_res );
void print_errors( GRID u, GRID f, GRID uex, int it );
void print_header( int lmax, int nc, int nu1, int nu2, int nuc, CYCLE cycle );
double residual_norm( GRID u, GRID f );
void residual( GRID u, GRID f, GRID d );
void interpolation( GRID uf, GRID uc );
void vcycle( GRIDLIST u, GRIDLIST f, GRIDLIST w, int l, int nu1, int nu2, int nuc );
void mgrid( GRIDLIST u, GRIDLIST f, GRIDLIST w, int l, int nu1, int nu2, int nuc, CYCLE cycle );

double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1.e-6;
}

void init_rhs( GRID u ) {
  //double x, y;
#pragma omp parallel for
  for (int j=0; j<=u.n; j++) {
    double y = j*u.h;
    for (int i=0; i<=u.n; i++) {
      double x=i*u.h;
      u.p[i][j] = 0.0;
      //u.p[i][j] = /*0.0;*/ -2.0*sin(x+y);
    }
  }
}    

void init_border( GRID u ) {
  double x, y;
  for (int i=0; i<=u.n+1; i++) {
    x = i*u.h;
    u.p[i][0] = 0.0; u.p[i][u.n+1] = 0.0;
    //u.p[i][0] = x*x; u.p[i][u.n+1] = 1.0+x*x;
    //u.p[i][0] = sin(x); u.p[i][u.n+1] = sin(1.0+x);
    //u.p[i][0] = x; u.p[i][u.n+1] = 1.0+x;
  }
  for (int j=1; j<=u.n; j++) {
    y = j*u.h;
    u.p[0][j] = 0.0; u.p[u.n+1][j] = 0.0;
    //u.p[0][j] = y*y; u.p[u.n+1][j] = y*y+1.0;
    //u.p[0][j] = sin(y); u.p[u.n+1][j] = sin(1.0+y);
    //u.p[0][j] = y; u.p[u.n+1][j] = 1.0+y;
  }
}

void init_solution( GRID u ) {
  //double x, y;
#pragma omp parallel for
  for (int j=0; j<=u.n+1; j++) {
    double y= j*u.h;
    for (int i=0; i<=u.n+1; i++) {
      double x=i*u.h;
      u.p[i][j] = 0.0;
      //u.p[i][j] = x*x+y*y;
      //u.p[i][j] = sin(x+y);
      //u.p[i][j] = x+y;
    }
  }
}

void create_grid( GRID *grid, int n ) {
  grid->n = n;
  grid->h = 1.0/(double)(n+1);
  grid->lda = n+2;//+(n+2)|64;

#ifdef USE_MM_MALLOC
  grid->p = (double **) _mm_malloc((n+2)*sizeof(double *),64);
  grid->p[0] = (double *) _mm_malloc((n+2)*(grid->lda)*sizeof(double),64);
#else
  grid->p = (double **) malloc((n+2)*sizeof(double *));
  grid->p[0] = (double *) malloc((n+2)*(grid->lda)*sizeof(double));
#endif
  for (int i=1; i<=n+1; i++)
    grid->p[i] = grid->p[i-1] + grid->lda;

  zero_grid(*grid);
}

void free_grid( GRID *grid ) {
#ifdef USE_MM_MALLOC
  _mm_free(grid->p[0]);
  _mm_free(grid->p);
#else
  free(grid->p[0]);
  free(grid->p);
#endif
}

void create_grid_list( GRIDLIST *list, int n, int lmax ) {
  *list = (GRID *) calloc(lmax+1,sizeof(GRID));
  if(*list == NULL) printf("create_grid_list: not enough memory");
  
  for (int i=0; i<=lmax; i++) {
    create_grid((*list)+i,n);
    n = 2*n+1;
  }
}

void free_grid_list( GRIDLIST *list, int lmax ) { 
  for (int i=0; i<=lmax; i++)
    free_grid((*list)+i);
  free(*list);
}

void print_grid( GRID u, char *c ) {
/*
  printf("\n GRID %s\n",c);
  for (int j=u.n+1; j>=0; j--) {
    printf("\n row %d:\n   ",j);
    for (int i=0; i<=u.n +1; i++)
      printf("%f ",u.p[i][j]);
  }
  printf("\n");
*/

    printf("\n-------------------------\n");
    for(int i = 0; i < 5; i++){
        for(int j = 0; j < 5; j++)
            printf("%6e\t", u.p[i][j]);
        printf("\t...,\t");
        for(int j = u.n+1-4; j <= u.n+1; j++)
            printf("%6e\t", u.p[i][j]);
        printf("\n");
    }
    printf("\t...,\n");
    for(int i = u.n+1-4; i <= u.n+1; i++){
        for(int j = 0; j < 5; j++)
            printf("%6e\t", u.p[i][j]);
        printf("\t...,\t");
        for(int j = u.n+1-4; j <= u.n+1; j++)
            printf("%6e\t", u.p[i][j]);
        printf("\n");
    }
}

void zero_grid( GRID u ) {
#pragma omp parallel for
  for (int i=0; i<= u.n+1; i++)
    for(int j=0; j<= u.n+1; j++) 
      u.p[i][j] = 0;
}

void random_grid( GRID u ) {
#pragma omp parallel for
  for (int i=0; i<= u.n+1; i++)
    for(int j=0; j<= u.n+1; j++) 
      u.p[i][j] = rand() / RAND_MAX;
}

void init_grid( GRID u, double a ) {
#pragma omp parallel for
  for (int i=0; i<=u.n+1; i++)
    for(int j=0; j<=u.n+1; j++) 
      u.p[i][j] = a;
}

double error( GRID u1, GRID u2 ) {
  double tmp = 0.0;
#pragma omp parallel for reduction(+:tmp) 
  for (int i=1; i<=u1.n; i++)
    for(int j=1; j<=u1.n; j++)
      tmp += (u1.p[i][j]-u2.p[i][j])*(u1.p[i][j]-u2.p[i][j]);
  tmp = sqrt(tmp)/u1.n;
  return tmp;
}

void correct( GRID u, GRID c ) {
#pragma omp parallel for
  for (int i=1; i<=u.n; i++)
    for (int j=1; j<=u.n; j++)
      u.p[i][j] += c.p[i][j];
}

void get_options( int *iters, int *nc, int *lmax, int *nu1, int *nu2, int *nuc, CYCLE *cycle, int *print_res ) {
  char cyc[2];
  char inp[2];
  int e;
  printf("# number of (internal) gridlines on coarsest level: ");
  e = scanf("%d",nc);
  printf("# additional finer levels: ");
  e = scanf("%d",lmax);
  printf("# multigrid cycle type (V), (W): " );
  e = scanf("%1s",cyc);
  switch(cyc[0]) {
  case 'V': case 'v': *cycle=V; break;
  case 'W': case 'w': *cycle=W; break;
  case 'S': case 's': *cycle=S; break;
  default  :  printf("wrong cycle type\n");
  }
  printf("# pre-smoothing steps:  ");
  e = scanf("%d",nu1);
  printf("# post-smoothing steps: ");
  e = scanf("%d",nu2);
  printf("# coarse-grid relaxation steps: ");
  e = scanf("%d",nuc);
  printf("# max iterations on finest level: ");
  e = scanf("%d",iters);
  printf("# time ? (y/n) : ");
  e = scanf("%1s", inp);
  switch(inp[0]) {
  case 'y': case 'Y': *print_res = 0;  break;
  case 'n': case 'N': *print_res = 1;  break;
  default : printf("wrong input");
  }
}

double oldresid,olderr;
void print_errors( GRID u, GRID f, GRID uex, int it ) {
  double rhoresid,rhoerr;
  double resid = residual_norm(u,f);
  double err   = error(u,uex);
  if (it==0) {
    printf("\n# %3dx%-3d :  %-3d :   %7e   %7e\n",u.n,u.n,it,err,resid);
  } else {
    rhoresid = (oldresid==0.0) ? 1.0 : resid/oldresid;
    rhoerr   = (olderr  ==0.0) ? 1.0 : err/olderr;
    printf("# %3dx%-3d :  %-3d :   %7e   %7e  :  %6f     %6f \n",u.n,u.n,it, err,resid,rhoerr,rhoresid);
  }
  oldresid = resid;
  olderr   = err;
}
  
void print_header( int lmax, int nc, int nu1, int nu2, int nuc, CYCLE cycle ) {
  char cyc[] = {'V', 'W', 'S'};
  printf(": levels: 0..%-2d, coarse grid: %3dx%-3d, %c cycle \n",lmax,nc,nc,cyc[(int)cycle]);
  printf("# pre-smoothing: %-2d, post-smoothing: %-2d, coarse relaxation: %-2d", nu1,nu2,nuc);
  printf("\n\n# discr     iter       error         residual       rho-error    rho-residual\n" );
  printf(   "# ------     ----       -----         --------       ---------    ------------\n" );
}

double residual_norm( GRID u, GRID f ) {
  double tmp, res = 0.0;  
  double invhh = 1.0/(u.h*u.h);
#pragma omp parallel for reduction(+:res)
  for (int i=1; i<=u.n; i++)
    for(int j=1; j<=u.n; j++) {
      tmp = (u.p[i-1][j]+u.p[i][j-1]+u.p[i+1][j]+u.p[i][j+1]-4.0*u.p[i][j])*invhh-f.p[i][j];
      res += tmp*tmp;
    }
  res = sqrt(res)/u.n;
  return res;
}

void residual( GRID u, GRID f, GRID d ) {
  double invhh = 1.0/(u.h*u.h);
#pragma omp parallel for collapse(2)
  for (int i=1; i<=u.n; i++)
    for (int j=1; j<=u.n; j++)
      d.p[i][j]=f.p[i][j]-(4.0*u.p[i][j]-u.p[i-1][j]-u.p[i+1][j]-u.p[i][j-1]-u.p[i][j+1])*invhh;
}

void restriction( GRID uf, GRID uc ) {
#pragma omp parallel for collapse(2)
  for (int i=2; i<=uf.n; i+=2)
    for (int j=2; j<=uf.n; j+=2)
      uc.p[i/2][j/2] = 0.0625*(uf.p[i-1][j-1]+uf.p[i+1][j-1]+uf.p[i-1][j+1]+uf.p[i+1][j+1])
	+0.125*(uf.p[i-1][j]+uf.p[i+1][j]+uf.p[i][j-1]+uf.p[i][j+1])
	+0.25*uf.p[i][j];
}

void interpolation( GRID uf, GRID uc ) {
#pragma omp parallel for collapse(2)
  for (int i=0; i<=uc.n; i++) {
    for (int j=0; j<=uc.n; j++) {
    /*
      if (i>0 && j>0) uf.p[2*i][2*j] = uc.p[i][j];
      if (j>0) uf.p[2*i+1][2*j] = uc.p[i][j];
      if (i>0) uf.p[2*i][2*j+1] = uc.p[i][j];
      uf.p[2*i+1][2*j+1] = uc.p[i][j];
    */
      if (i>0 && j>0) uf.p[2*i][2*j] = uc.p[i][j];
      if (j>0) uf.p[2*i+1][2*j] = (uc.p[i+1][j] + uc.p[i][j]) / 2.0;
      if (i>0) uf.p[2*i][2*j+1] = (uc.p[i][j+1] + uc.p[i][j]) / 2.0;
      uf.p[2*i+1][2*j+1] = (uc.p[i][j] + uc.p[i][j+1] + uc.p[i+1][j] + uc.p[i+1][j+1]) / 4.0;
    }
  }
}

/*
void vcycle( GRIDLIST u, GRIDLIST f, GRIDLIST w, int l, int nu1, int nu2, int nuc ) {
  if (l==0) {
    //u[0].p[1][1] = (u[0].h)*(u[0].h) / 4.0 * f[0].p[1][1];
    jacobi(&u[0],&f[0],&w[0],nuc);
  } else {
    jacobi(&u[l],&f[l],&w[l],nu1);
    residual(u[l],f[l],w[l]);    
    restriction(w[l],f[l-1]);
    zero_grid(u[l-1]);
    vcycle(u,f,w,l-1,nu1,nu2,nuc);
    interpolation(w[l],u[l-1]);
    correct(u[l],w[l]);
    jacobi(&u[l],&f[l],&w[l],nu2);
  }
}
*/

void mgrid( GRIDLIST u, GRIDLIST f, GRIDLIST w, int l, int nu1, int nu2, int nuc, CYCLE cycle) {
  if (l==0) {
    //u[0].p[1][1][1] = (u[0].h)*(u[0].h) / 6.0 * f[0].p[1][1][1];
    jacobi(&u[0],&f[0],&w[0],nuc);
  } else {
    jacobi(&u[l],&f[l],&w[l],nu1);
    residual(u[l],f[l],w[l]);    
    restriction(w[l],f[l-1]);

    if (cycle == V) {
      zero_grid(u[l-1]);
      mgrid(u,f,w,l-1,nu1,nu2,nuc,cycle);
    } else if (cycle == W) {
      zero_grid(u[l-1]);
      mgrid(u,f,w,l-1,nu1,nu2,nuc,cycle);
      mgrid(u,f,w,l-1,nu1,nu2,nuc,cycle);
    } else if (cycle == S) {
      if (l==1) {
	zero_grid(u[l-1]);
	mgrid(u,f,w,l-1,nu1,nu2,nuc,cycle);
      } else {
	restriction(f[l-1],f[l-2]);
	zero_grid(u[l-2]);
	mgrid(u,f,w,l-2,nu1,nu2,nuc,cycle);
	interpolation(u[l-1],u[l-2]);
      }
    }

    interpolation(w[l],u[l-1]);
    correct(u[l],w[l]);
    jacobi(&u[l],&f[l],&w[l],nu2);
  }
}

int main( void ) {
  int lmax,nc,nu1,nu2,nuc,iters,print_res;
  double tstart, ttotal,tmin=10000000.0;
  CYCLE cycle;
  GRIDLIST u,f,w,uex;
    
  get_options(&iters,&nc,&lmax,&nu1,&nu2,&nuc,&cycle,&print_res);
  print_header(lmax,nc,nu1,nu2,nuc,cycle);
  
  create_grid_list(&uex,nc,lmax);
  create_grid_list(&u,nc,lmax);
  create_grid_list(&f,nc,lmax);
  create_grid_list(&w,nc,lmax);

  ttotal = 0;
  int s;
  for(s=0;s<1;s++){
  init_grid(u[lmax], 1.0);
  init_border(u[lmax]);
  init_border(w[lmax]);
  init_rhs(f[lmax]);
  init_solution(uex[lmax]);
  print_errors(u[lmax],f[lmax],uex[lmax],0);
  tstart = get_time();
  int i = 0;
    char ch = 'u';

  //print_grid(u[lmax], &ch);
  while (i < iters && error(u[lmax],uex[lmax]) > 1e-12) {  
    mgrid(u,f,w,lmax,nu1,nu2,nuc,cycle);
    //vcycle(u,f,w,lmax,nu1,nu2,nuc);

    //print_grid(u[lmax], &ch);
    i++;
    if (print_res) 
      print_errors(u[lmax],f[lmax],uex[lmax],i);
  }
  ttotal += (get_time()-tstart);
  if(ttotal<tmin) tmin = ttotal;
 }
    //print_grid(u[lmax], &ch);
  if (!print_res)
    //printf("\n# total time: \n%7f\t%i #sec.\n", ttotal, i);
    printf("\n%7f ms\n", tmin*1000);
  free_grid_list(&uex,lmax);
  free_grid_list(&u,lmax);
  free_grid_list(&f,lmax);
  free_grid_list(&w,lmax);
}


