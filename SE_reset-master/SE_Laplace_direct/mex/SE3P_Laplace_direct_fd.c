/*
Equivalent MATLAB code (for potential):

function u = SE3P_Laplace_direct_fd(idx, x, f, opt)
N = numel(f);
xi = opt.xi;
xi2 = xi*xi;
TwoPiOverL = 2*pi./opt.box;
c = 4*pi/(opt.box(1)*opt.box(2)*opt.box(3));
a = 1/(4*xi2);
u = zeros(N,1);
for target=idx
  xt = x(target,:);
  u_t = 0;
  for j1=-opt.layers:opt.layers
    for j2=-opt.layers:opt.layers
      for j3=-opt.layers:opt.layers
        if j1 == 0 && j2 == 0 && j3 == 0
          continue
        end
        k = TwoPiOverL .* [j1 j2 j3];
        k2 = dot(k,k);
        z = 0;
        for source=1:N
          X = xt - x(source,:);
          fn = f(source);
          z = z + fn*cos(-dot(k,X)); % real part of exp(-i*dot(k,X))
                                     % (imaginary part cancels)
        end
        u_t = u_t + c * z * exp(-a*k2)/k2;
      end
    end
  end
  u(target) = u_t;
end
*/

#include "mex.h"
#include "SE_direct.h"

#define IDX prhs[0]
#define X   prhs[1] // Source locations
#define Q   prhs[2] // Source strengths
#define OPT prhs[3] // Parameters

#define OUT plhs[0] // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

/* common option-unpacking */
void unpack_opt(ewald_opts* opt, const mxArray* mx_opt)
{
    // mandatory options -- will trigger core dump if missing
    opt->xi = mxGetScalar(mxGetField(mx_opt,0,"xi"));
    double* box =  mxGetPr(mxGetField(mx_opt,0,"box"));

    opt->box[0] = box[0];
    opt->box[1] = box[1];
    opt->box[2] = box[2];

    // layers: mandatory for ewald sums that are truncated 
    const mxArray* mx_layers = mxGetField(mx_opt,0,"layers");
    if (mx_layers)
        opt->layers = (int)mxGetScalar(mx_layers);
    else
        opt->layers = -1;

    // rc: mandatory for short-range real sum 
    const mxArray* mx_rc = mxGetField(mx_opt,0,"real_cutoff");
    if (mx_rc)
        opt->rc = mxGetScalar(mx_rc);
    else
        opt->rc = -1;
}

// MATLAB (one-based, doubles) to C (zero-based, integers) index translation
void index_translation(int* idx, const double* idx_d, int N)
{
    for (int i=0; i<N; i++)
        idx[i] = (int)idx_d[i] - 1;
}


/* no input checking is done */
void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
    // input dims
    const int N = mxGetM(X);

    const int num_eval = mxGetNumberOfElements(IDX);
    const double* idx_d = mxGetPr(IDX);
    int* idx = mxMalloc(num_eval*sizeof(int));
    index_translation(idx, idx_d, num_eval);

    const double* x = mxGetPr(X);
    const double* q = mxGetPr(Q);

#ifndef FORCE
    // Allocate matrix for the potential
    OUT = mxCreateDoubleMatrix(num_eval, 1, mxREAL);
    double* restrict pot = mxGetPr(OUT);
#else
    // Allocate matrix for the force vectors
    OUT = mxCreateDoubleMatrix(num_eval, 3, mxREAL);
    double* restrict force = mxGetPr(OUT);
#endif

    ewald_opts opt;
    unpack_opt(&opt, OPT);

    if (VERBOSE)
    {
        mexPrintf("[EWALD (%s)] MEX N=(%d,%d) ","FS3P", N, num_eval);
        mexPrintf("xi = %.2f [rc = %.2f, layers=%d]\n",
                  opt.xi, opt.rc, opt.layers);
    }

    // main computation
    double k[3];
    double k2, z;
    double c = 4*PI/(opt.box[0]*opt.box[1]*opt.box[2]);
    double fac[3] = {2.*PI/opt.box[0], 2.*PI/opt.box[1], 2.*PI/opt.box[2]};
    double a = 1.0/(4*opt.xi*opt.xi);
#ifdef _OPENMP
#pragma omp parallel for private(k, k2, z)
#endif
    for (int m=0; m<num_eval; m++)
    {
#ifndef FORCE
        double p = 0;
#else
        double f[] = {0, 0, 0};
#endif
        double xm[3] = {x[idx[m]],x[idx[m]+N],x[idx[m]+2*N]};
        for (int j0 = -opt.layers; j0<=opt.layers; j0++)
            for (int j1 = -opt.layers; j1<=opt.layers; j1++)
                for (int j2 = -opt.layers; j2<=opt.layers; j2++)
                {
                    if (j0 == 0 && j1 == 0 && j2==0)
                        continue;

                    k[0] = fac[0]*j0;
                    k[1] = fac[1]*j1;
                    k[2] = fac[2]*j2;
                    k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];

                    z=0;
                    double tmp;
                    for (int n = 0; n<N; n++)
                    {
                        tmp = -(k[0]*(xm[0]-x[n]    )+
                                k[1]*(xm[1]-x[n+N]  )+
                                k[2]*(xm[2]-x[n+2*N])
                                );
/* Only the real part remains here, since the imaginary part cancels */
#ifndef FORCE
                        z += q[n]*cos(tmp);
#else
                        z += q[n]*sin(tmp);
#endif
                    }
#ifndef FORCE
                    p += z*exp(-a*k2)/k2;
#else
                    tmp = z*exp(-a*k2)/k2;
                    f[0] += k[0]*tmp;
                    f[1] += k[1]*tmp;
                    f[2] += k[2]*tmp;
#endif
                }
#ifndef FORCE
        pot[m] += c*p;
#else
        force[m           ] += c*f[0];
        force[m+  num_eval] += c*f[1];
        force[m+2*num_eval] += c*f[2];
#endif
    }
    mxFree(idx);
}
