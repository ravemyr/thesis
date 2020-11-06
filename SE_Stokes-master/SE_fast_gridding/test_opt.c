#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "test_opt.h"


/* This routine can be used for optimiziation purpose  */
/* on kernels in SE_fgg.c. */

int main(){
  const int N=100000;
  const int P = 24;
  const double L=1;
  const double xi=4;
  const double beta = 6.*P/2;
  int M[3] = {28,28,28};
  double box[3]={L,L,L};
  double h = L/(double) M[0];
  double w = P*h/2.;
  double m = .98*sqrt(PI*P);
  double eta = pow((2.*w*xi/m),2.);
  double c = 2.*xi*xi/eta;
  
  // allocate memory
  double * restrict x = (double*) malloc(sizeof(double)*N*3);
  double * restrict q = (double*) malloc(sizeof(double)*N);
  double * restrict p = (double*) malloc(sizeof(double)*N);

  // fill arrays
  double average = 0;
  srand(0);
  for (int b=0; b<N; b++) {
    for (int d=0; d<3; d++) {
      x[d*N+b] = rand()/(double) RAND_MAX * L;
    }
    q[b] = rand()/(double) RAND_MAX - .5;
    average += q[b];
  }
  average /= N;
  for (int b=0; b<N; b++) {
    q[b] -= average;
  }
  
  // pack parameters
  SE_FGG_params params;
  set_params(&params, N, box, M, P, c, beta);
  
  // scratch arrays
  SE_FGG_work work;
  SE_FGG_allocate_workspace(&work, &params,1,0);
  
  // coordinates and charges
  SE_state st = {.x = x,  .q = q};
    
  // now do the gridding
  SE_FGG_grid_kaiser(&work, &st, &params);

  // now do the gathering
  st.q = NULL;
  SE_FGG_int_kaiser(p, &work, &st, &params);

  printf("Finished !  \n");
  
  return 0;
}
