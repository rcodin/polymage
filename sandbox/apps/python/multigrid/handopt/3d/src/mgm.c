#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "grid.h"
#include "intergrid.h"
#include "smoothers.h"

typedef enum {V,W,S} CYCLE;
typedef enum {JACOBI,CHEBYSHEV} METHOD;

void get_options( int *iters, int *nc, double *mineigen, int *lmax,
		  int *nu1, int *nu2, int *nuc, CYCLE *cycle, METHOD *smoother, int *print_res );
double residual( GRID u, GRID f );
void print_errors( GRID u, GRID f, GRID uex, int it );
void print_header( double mineigen, int lmax, int nc, int nu1, int nu2, int nuc, CYCLE cycle, METHOD smoother );
long long int defect( GRID u, GRID f, GRID d );
long long int smooth( METHOD smoother, GRID &u, GRID &f, GRID &w, int nu, double mineigen );
long long int mgrid( GRIDLIST u, GRIDLIST f, GRIDLIST w, int l, int nu1, int nu2, int nuc, CYCLE cycle, METHOD smoother, double mineigen );
/*
void defect( GRID u, GRID f, GRID d );
void smooth( METHOD smoother, GRID &u, GRID &f, GRID &w, int nu, double mineigen );
void mgrid( GRIDLIST u, GRIDLIST f, GRIDLIST w, int l, int nu1, int nu2, int nuc, CYCLE cycle, METHOD smoother, double mineigen );
*/

void get_options( int *iters, int *nc, double *mineigen, int *lmax,
		  int *nu1, int *nu2, int *nuc, CYCLE *cycle, METHOD *smoother, int *print_res ) {
  
  char cyc[2],inp[2];

  printf("\n\n# ITERATIEVE METHODES voor het OPLOSSEN van de POISSON VERGELIJKING\n\n");
  printf("\n# aantal (inwendige) roosterlijnen op grofste rooster: ");
  scanf("%d",nc);
  printf("# aantal bijkomende fijne roosters: ");
  scanf("%d",lmax);
  printf("# multigrid cycle type (V), (W): " );
  scanf("%1s",cyc);
  switch(cyc[0]) {
  case 'V': case 'v': *cycle=V; break;
  case 'W': case 'w': *cycle=W; break;
  case 'S': case 's': *cycle=S; break;
  default  :  print_msg("verkeerd cycle type");
  }
  printf("# multigrid smoother (J for Jacobi), (C for Chebyshev): " );
  scanf("%1s",cyc);
  switch(cyc[0]) {
  case 'J' : case 'j' :  *smoother=JACOBI; break;
  /* case 'C' : case 'c' :  { */
  /*   *smoother=CHEBYSHEV;  */
  /*   printf("# l_min = c l_max, geef c:  "); */
  /*   scanf("%lf",mineigen);     */
  /*   break; */
  /* } */
  default  :  print_msg("verkeerde smoother type");
  }
  printf("# aantal pre-smoothing stappen:  ");
  scanf("%d",nu1);
  printf("# aantal post-smoothing stappen: ");
  scanf("%d",nu2);
  printf("# aantal grof-rooster-relaxatie stappen: ");
  scanf("%d",nuc);
  printf("# aantal iteraties op het fijnste rooster: ");
  scanf("%d",iters);
  
  printf("# Tijden opmeten ? (j/n) : ");
  scanf("%1s", inp);
  switch(inp[0]) {
  case 'j': case 'J': *print_res = 0;  break;
  case 'n': case 'N': *print_res = 1;  break;
  default : print_msg("verkeerde invoer");
  }
}

double residual( GRID u, GRID f ) {
  if (u.n != f.n) print_msg("residual: wrong grid-sizes");
  double tmp, res = 0.0;  
  double invhh = 1.0/(u.h*u.h);
#pragma omp parallel for collapse(3) reduction(+:res) private(tmp)
  for (int i=1; i<=u.n; i++)
    for (int j=1; j<=u.n; j++)
      for (int k=1; k<=u.n; k++) {
	tmp = (u.p[i-1][j][k]+
	       u.p[i+1][j][k]+	       
	       u.p[i][j-1][k]+	       
	       u.p[i][j+1][k]+
	       u.p[i][j][k-1]+
	       u.p[i][j][k+1]-
	       6.0*u.p[i][j][k])*invhh-f.p[i][j][k];
	res += tmp*tmp;
      }
  res = sqrt(res)/u.n;
  return res;
}

double oldresid,olderr;
void print_errors( GRID u, GRID f, GRID uex, int it ) {
  double rhoresid,rhoerr;
  double resid = residual(u,f);
  double err   = error(u,uex);

  if (it==0) {
    printf("\n# %dx%dx%d :  %-3d :   %20e   %20e\n",u.n,u.n,u.n,it,err,resid);
  } else {
    rhoresid = (oldresid==0.0) ? 1.0 : resid/oldresid;
    rhoerr   = (olderr  ==0.0) ? 1.0 : err/olderr;
    printf("# %dx%dx%d :  %-3d :   %20e   %20e  :  %6f     %6f \n",u.n,u.n,u.n,it, err,resid,rhoerr,rhoresid);
  }
  oldresid = resid;
  olderr   = err;
}
  
void print_header( double mineigen, int lmax, int nc, int nu1, int nu2, int nuc, CYCLE cycle, METHOD smoother ) {
  char cyc[] = {'V','W','S'};
  char* str[] = {"JACOBI", "CHEBYSHEV"};
  printf("\n\n\n# %s smoother ", str[(int)smoother]);
  printf(": levels: 0..%-2d, coarse grid: %dx%dx%d, %c cycle, lambda_min= %f lambda_max\n",lmax,nc,nc,nc,cyc[(int)cycle],mineigen);
  printf("# pre-smoothing: %-2d, post-smoothing: %-2d, coarse relaxation: %-2d", nu1,nu2,nuc);
  printf("\n\n# discr     iter       error         residual       rho-error    rho-residual\n" );
  printf(   "# ------     ----       -----         --------       ---------    ------------\n" );
}

//void defect( GRID u, GRID f, GRID d ) {
long long int defect( GRID u, GRID f, GRID d ) {
    long long int count = 0;
  if (u.n != f.n || u.n != d.n) print_msg("residual: wrong grid sizes");
  double invhh = 1.0/(u.h*u.h);
#pragma omp parallel for collapse(3)
  for (int i=1; i<=u.n; i++)
    for (int j=1; j<=u.n; j++)
      for (int k=1; k<=u.n; k++)
    
	d.p[i][j][k]=f.p[i][j][k]-(6.0*u.p[i][j][k]-
				   u.p[i-1][j][k]-
				   u.p[i+1][j][k]-
				   u.p[i][j-1][k]-
				   u.p[i][j+1][k]-
				   u.p[i][j][k-1]-
				   u.p[i][j][k+1] )*invhh;
    

   // count++;

    return count;
}

//void smooth( METHOD smoother, GRID &u, GRID &f, GRID &w, int nu, double mineigen ) {
long long int smooth( METHOD smoother, GRID &u, GRID &f, GRID &w, int nu, double mineigen ) {
  double lmax = 8.0*(u.n+1)*(u.n+1);  
  double lmin = mineigen * lmax;

    long long int count = 0;
  //if (smoother == JACOBI)
  count = jacobi(u,f,w,nu);
  //else if (smoother == CHEBYSHEV)
  //chebyshev(u,f,w,nu,lmin,lmax);

    return count;
}

//void mgrid( GRIDLIST u, GRIDLIST f, GRIDLIST w, int l, int nu1, int nu2, int nuc, CYCLE cycle, METHOD smoother, double mineigen ) {
long long int mgrid( GRIDLIST u, GRIDLIST f, GRIDLIST w, int l, int nu1, int nu2, int nuc, CYCLE cycle, METHOD smoother, double mineigen ) {

    long long int count = 0;
  if (l==0) {
    u[0].p[1][1][1] = (u[0].h)*(u[0].h) / 6.0 * f[0].p[1][1][1];
    count += smooth(smoother,u[0],f[0],w[0],nuc,mineigen);
  } else {
    count += smooth(smoother,u[l],f[l],w[l],nu1,mineigen);
    count += defect(u[l],f[l],w[l]);    
    count += restriction(w[l],f[l-1]);

    if (cycle == V) {
      zero_grid(u[l-1]);
      count += mgrid(u,f,w,l-1,nu1,nu2,nuc,cycle,smoother,mineigen);
    } else if (cycle == W) {
      zero_grid(u[l-1]);
      count += mgrid(u,f,w,l-1,nu1,nu2,nuc,cycle,smoother,mineigen);
      count += mgrid(u,f,w,l-1,nu1,nu2,nuc,cycle,smoother,mineigen);
    } else if (cycle == S) {
      if (l==1) {
	zero_grid(u[l-1]);
	mgrid(u,f,w,l-1,nu1,nu2,nuc,cycle,smoother,mineigen);
      } else {
	count += restriction(f[l-1],f[l-2]);
    zero_grid(u[l-2]);
	count += mgrid(u,f,w,l-2,nu1,nu2,nuc,cycle,smoother,mineigen);
	count += prolong(u[l-1],u[l-2]);
      }
    }

    count += prolong(w[l],u[l-1]);
    count += correct(u[l],w[l]);
    count += smooth(smoother,u[l],f[l],w[l],nu2,mineigen);
  }

    return count;
}

int main( void ) {
  int lmax,nc,nu1,nu2,nuc,iters,print_res;
  double tstart, ttotal, tmin;
  double mineigen;  
  CYCLE cycle;
  METHOD smoother;
  GRIDLIST u,f,w,uex;
    
  get_options(&iters,&nc,&mineigen,&lmax,&nu1,&nu2,&nuc,&cycle,&smoother,&print_res);
  print_header(mineigen,lmax,nc,nu1,nu2,nuc,cycle,smoother);
  
  create_grid_list(&uex,nc,lmax);
  create_grid_list(&u,nc,lmax);
  create_grid_list(&f,nc,lmax);
  create_grid_list(&w,nc,lmax);

    long long int count = 0;

  tmin = 1000.0;
  int s, i;
  for (s=0; s<1; s++){
    init_grid(u[lmax], 1.0);
    init_border(u[lmax]);
    init_border(w[lmax]);
    init_rhs(f[lmax]);
    init_solution(uex[lmax]);
    print_errors(u[lmax],f[lmax],uex[lmax],0);
    tstart = get_time();
    i=0;
    while (i < iters && error(u[lmax],uex[lmax]) > 1e-12) {  
      count += mgrid(u,f,w,lmax,nu1,nu2,nuc,cycle,smoother,mineigen);
      i++;
      if (print_res) 
	print_errors(u[lmax],f[lmax],uex[lmax],i);
    }
    ttotal = (get_time()-tstart);
    if (ttotal < tmin) tmin = ttotal;
  }
  if (!print_res)
    printf("\n# totale berekeningstijd: \n%7f\t%i #sec.\n\n", ttotal, i);
  free_grid_list(&uex,lmax);
  free_grid_list(&u,lmax);
  free_grid_list(&f,lmax);
  free_grid_list(&w,lmax);
    printf("%lld\n", count);
}


