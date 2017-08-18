#include "intergrid.h"

//void restriction( GRID uf, GRID uc ) {
long long int restriction( GRID uf, GRID uc ) {
    long long int count = 0;
  if (uf.n != 2*uc.n+1) print_msg("restriction: wrong grid sizes");
#pragma omp parallel for collapse(2)
  for (int i=2; i<=uf.n; i+=2)
    for (int j=2; j<=uf.n; j+=2)
      for (int k=2; k<=uf.n; k+=2)
   
	uc.p[i/2][j/2][k/2] = 1.0/64.0* (uf.p[i-1][j-1][k-1] + uf.p[i+1][j-1][k-1] + uf.p[i-1][j+1][k-1] + uf.p[i+1][j+1][k-1] +
				    uf.p[i-1][j-1][k+1] + uf.p[i+1][j-1][k+1] + uf.p[i-1][j+1][k+1] + uf.p[i+1][j+1][k+1])
	  +1.0/32.0* (uf.p[i-1][j][k-1] + uf.p[i+1][j][k-1] + uf.p[i][j-1][k-1] + uf.p[i][j+1][k-1] +
		      uf.p[i-1][j-1][k] + uf.p[i-1][j+1][k] + uf.p[i+1][j+1][k] + uf.p[i+1][j-1][k] +
		      uf.p[i-1][j][k+1] + uf.p[i+1][j][k+1] + uf.p[i][j-1][k+1] + uf.p[i][j+1][k+1])
	  +1.0/16.0* (uf.p[i-1][j][k] + uf.p[i+1][j][k] + uf.p[i][j-1][k] + uf.p[i][j+1][k] + uf.p[i][j][k-1] + uf.p[i][j][k+1])
	  +1.0/8.0* uf.p[i][j][k];
    
    //count++;

    return count;
}

//void prolong( GRID uf, GRID uc ) {
long long int prolong( GRID uf, GRID uc ) {
    long long int count = 0;
  if (uf.n != 2*uc.n+1) print_msg("prolong: wrong grid sizes");

#pragma omp parallel for collapse(2)
  for (int i=0; i<=uc.n; i++) {
    for (int j=0; j<=uc.n; j++) {
      for (int k=0; k<=uc.n; k++) {
    
	if (i>0 && j>0 && k>0) uf.p[2*i][2*j][2*k] = uc.p[i][j][k];
	if (j>0 && k>0) uf.p[2*i+1][2*j][2*k] = (uc.p[i+1][j][k] + uc.p[i][j][k]) / 2.0;
	if (i>0 && k>0) uf.p[2*i][2*j+1][2*k] = (uc.p[i][j+1][k] + uc.p[i][j][k]) / 2.0;
	if (i>0 && j>0) uf.p[2*i][2*j][2*k+1] = (uc.p[i][j][k+1] + uc.p[i][j][k]) / 2.0;

	if (i>0) uf.p[2*i][2*j+1][2*k+1] = (uc.p[i][j][k] + uc.p[i][j+1][k] + uc.p[i][j][k+1] + uc.p[i][j+1][k+1]) / 4.0;
	if (j>0) uf.p[2*i+1][2*j][2*k+1] = (uc.p[i][j][k] + uc.p[i+1][j][k] + uc.p[i][j][k+1] + uc.p[i+1][j][k+1]) / 4.0;
	if (k>0) uf.p[2*i+1][2*j+1][2*k] = (uc.p[i][j][k] + uc.p[i+1][j][k] + uc.p[i][j+1][k] + uc.p[i+1][j+1][k]) / 4.0;

	uf.p[2*i+1][2*j+1][2*k+1] = (uc.p[i][j][k] + uc.p[i][j+1][k] + uc.p[i+1][j][k] + uc.p[i+1][j+1][k] +
				     uc.p[i][j][k+1] + uc.p[i][j+1][k+1] + uc.p[i+1][j][k+1] + uc.p[i+1][j+1][k+1]) / 8.0;
    
   // count++;

      }
    }
  }
    return count;
}
