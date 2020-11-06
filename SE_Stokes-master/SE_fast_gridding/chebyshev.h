void interpolant ( double a, double b, int n, double c[],
		   int p, int p_half, double x , double *z);
double *chebyshev_coefficients ( double a, double b, int n, double ow2, double beta,
				 double f ( double x, double ow2, double beta ) );
double *chebyshev_interpolant ( double a, double b, int n, double c[], int m, 
  double x[] );
double *chebyshev_zeros ( int n );
void timestamp ( void );
