#ifndef grid_h_
#define grid_h_

//#include "utils.h"

typedef struct{ double ***p; int n; double h; int lda; } GRID;
typedef GRID *GRIDLIST;

void init_rhs( GRID u);
void init_border( GRID u );
void init_solution( GRID u );

void create_grid( GRID *grid, int n );
void free_grid( GRID *grid );
void create_grid_list( GRIDLIST *list, int n, int lmax );
void free_grid_list( GRIDLIST *list, int lmax );
void print_grid( GRID u );
void zero_grid( GRID u );
void init_grid( GRID u, double a );
void random_grid( GRID u );

double error( GRID u1, GRID u2 );
//void correct( GRID u, GRID c );
long long int correct( GRID u, GRID c );

#endif // grid_h_
