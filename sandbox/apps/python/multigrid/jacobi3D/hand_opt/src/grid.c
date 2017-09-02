#include "grid.h"
#include "utils.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

void init_rhs( GRID u ) {
#pragma omp parallel for
  for (int i=1; i<=u.n; i++)
    for (int j=1; j<=u.n; j++)
      for (int k=1; k<=u.n; k++)
	u.p[i][j][k] = 0.0;
}    

void init_border( GRID u ) {
  for (int i=0; i<=u.n+1; i++)
    for (int j=0; j<=u.n+1; j++) {
      u.p[i][j][0] = 0.0;
      u.p[i][j][u.n+1] = 0.0;
    }
  for (int i=0; i<=u.n+1; i++)
    for (int k=0; k<=u.n+1; k++) {
      u.p[i][0][k] = 0.0;
      u.p[i][u.n+1][k] = 0.0;
    }
  for (int j=0; j<=u.n+1; j++)
    for (int k=0; k<=u.n+1; k++) {
      u.p[0][j][k] = 0.0;
      u.p[u.n+1][j][k] = 0.0;
    }
}

void init_solution( GRID u ) {
  zero_grid(u);
}

void create_grid( GRID *grid, int n ) {
  grid->n = n;
  grid->h = 1.0/(double)(n+1);
  grid->lda = n+2;

#ifdef USE_MM_MALLOC
  grid->p = (double ***) _mm_malloc((n+2)*sizeof(double **),64);
  grid->p[0] = (double **) _mm_malloc((n+2)*(grid->lda)*sizeof(double *),64);
  grid->p[0][0] = (double *) _mm_malloc((n+2)*(grid->lda)*(grid->lda)*sizeof(double),64);
#else
  grid->p = (double ***) malloc((n+2)*sizeof(double **));
  grid->p[0] = (double **) malloc((n+2)*(grid->lda)*sizeof(double *));
  grid->p[0][0] = (double *) malloc((n+2)*(grid->lda)*(grid->lda)*sizeof(double));
#endif
  for (int i=0; i<=n+1; i++) {
    grid->p[i] = grid->p[0] + i*grid->lda;
    for (int j=0; j<=n+1; j++) {
      grid->p[i][j] = grid->p[0][0] + j*grid->lda + i*grid->lda*grid->lda;
    }
  }
  
  zero_grid(*grid);
}

void free_grid( GRID *grid ) {
#ifdef USE_MM_MALLOC
  _mm_free(grid->p[0][0]);
  _mm_free(grid->p[0]);
  _mm_free(grid->p);
#else
  free(grid->p[0][0]);
  free(grid->p[0]);
  free(grid->p);
#endif
}

void create_grid_list( GRIDLIST *list, int n, int lmax ) {
  *list = (GRID *) calloc(lmax+1,sizeof(GRID));
  if(*list == NULL) print_msg("create_grid_list: not enough memory");
  
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


void print_grid( GRID u ) {
  cout << "\n GRID \n";
  for (int i=0; i<=u.n+1; i++) {
    for (int j=0; j<=u.n+1; j++) {
      for (int k=0; k<=u.n+1; k++)
	cout << u.p[i][j][k] << "  ";
      cout << "\n";
    }
    cout << "\n";
  }
  cout << "\n";
}

/* void print_grid( GRID u, char *c ) { */
/*   printf("\n GRID %s\n",c); */
/*   for (int j=u.n+1; j>=0; j--) { */
/*     printf("\n row %d:\n   ",j); */
/*     for (int i=0; i<=u.n +1; i++) */
/*       printf("%f ",u.p[i][j]); */
/*   } */
/*   printf("\n"); */
/* } */

void zero_grid( GRID u ) {
#pragma omp parallel for
  for (int i=0; i<=u.n+1; i++)
    for (int j=0; j<=u.n+1; j++)
      for (int k=0; k<=u.n+1; k++) 
	u.p[i][j][k] = 0;
}

void random_grid( GRID u ) {
#pragma omp parallel for
  for (int i=0; i<=u.n+1; i++)
    for(int j=0; j<=u.n+1; j++) 
      for(int k=0; k<=u.n+1; k++) 
	u.p[i][j][k] = rand() / RAND_MAX;
}

void init_grid( GRID u, double a ) {
#pragma omp parallel for
  for (int i=0; i<=u.n+1; i++)
    for(int j=0; j<=u.n+1; j++) 
      for(int k=0; k<=u.n+1; k++) 
	u.p[i][j][k] = a;
}

double error( GRID u1, GRID u2 ) {
  if (u1.n != u2.n) print_msg("error: wrong grid-sizes");
  double tmp = 0.0;
#pragma omp parallel for reduction(+:tmp) 
  for (int i=1; i<=u1.n; i++)
    for(int j=1; j<=u1.n; j++)
      for(int k=1; k<=u1.n; k++)
	tmp += (u1.p[i][j][k]-u2.p[i][j][k])*(u1.p[i][j][k]-u2.p[i][j][k]);
  tmp = sqrt(tmp)/u1.n;
  return tmp;
}

//void correct( GRID u, GRID c ) {
long long int correct( GRID u, GRID c ) {
    long long int count = 0;
  if (u.n != c.n) print_msg("correct: wrong grid sizes");
#pragma omp parallel for
  for (int i=1; i<=u.n; i++)
    for (int j=1; j<=u.n; j++)
      for (int k=1; k<=u.n; k++)
	u.p[i][j][k] += c.p[i][j][k];
    //count++;

    return count;
}
