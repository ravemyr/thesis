#include <stdio.h>
#include <string.h>
#include <mex.h>
#include "SE_fgg.h"
#include "SE_fg_windows.c"


// Get field s from matlab struct p. abort if field is missing
void* get_arg(const mxArray* p, const char* s)
{
  mxArray* x=mxGetField(p,0,s);
  if (x) return mxGetData(x);
  else
  {
    mexErrMsgTxt("Missing mandatory parameter in struct");
    return (void*) NULL;
  }
}
void* get_arg_str(const mxArray* p, const char* s)
{
  mxArray* x=mxGetField(p,0,s);
  if (x) return mxArrayToString(x);
  else
  {
    mexErrMsgTxt("Missing mandatory parameter in struct");
    return (void*) NULL;
  }
}

KaiserlikeWindow get_kaiserlike_window(const char* window_id)
{
  if (strcmp(window_id, "expsemicirc") == 0) {
    return &SE_fg_window_expsemicirc;
  } else if (strcmp(window_id, "kaiser_exact") == 0) {
    return &SE_fg_window_kaiser_exact;
  } else if (strcmp(window_id, "kaiser_poly") == 0) {
    return &SE_fg_window_kaiser_poly;
  } else {
    mexErrMsgTxt("Unsupported window function");
    return (KaiserlikeWindow) NULL;
  }
}

KaiserlikeWindow get_kaiserlike_window_deriv(const char* window_id)
{
  if (strcmp(window_id, "kaiser_exact") == 0) {
    return &SE_fg_window_kaiser_exact_deriv;
  } else if (strcmp(window_id, "kaiser_poly") == 0) {
    return &SE_fg_window_kaiser_poly_deriv;
  } else {
    mexErrMsgTxt("Unsupported window function");
    return (KaiserlikeWindow) NULL;
  }
}

// get relevant parameters from mxArray and populate FGG params struct
void SE_FGG_MEX_params(SE_FGG_params* params, const mxArray* OPT, int N)
{

#ifdef THREE_PERIODIC

    const double* m    = (double*) get_arg(OPT,"M");
    const double* p    = (double*) get_arg(OPT,"P");
    const double* box  = (double*) get_arg(OPT,"box");
#ifndef KAISER
    const double* c    = (double*) get_arg(OPT,"c");
    params->c = *c;
    params->d = pow(params->c/PI,1.5);
#else
    const double* beta = (double*) get_arg(OPT,"beta");
    params->beta = *beta;
    char* window = (char*) get_arg_str(OPT,"window");
    params->window = get_kaiserlike_window(window);
#ifdef FORCE
    params->window_deriv = get_kaiserlike_window_deriv(window);
#endif
    if (strcmp(window, "kaiser_poly") == 0) {
      params->use_polynomial_window = 1;
      params->polynomial_degree = (int) * (double*) get_arg(OPT,"polynomial_degree");
    } else {
      params->use_polynomial_window = 0;
    }
    mxFree(window);
#endif

    params->N = N;
    params->P = (int) *p;
    params->P_half=half( (int) *p );
    params->P_eval_window = params->P;
    params->h = box[0]/m[0];
    params->a = -FGG_INF;

    params->dims[0] = (int) m[0];
    params->dims[1] = (int) m[1];
    params->dims[2] = (int) m[2];

    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1]+params->P;
    params->npdims[2] = params->dims[2]+params->P;

#endif

#ifdef TWO_PERIODIC

    const double* m    = (double*) get_arg(OPT,"M");
    const double* mz   = (double*) get_arg(OPT,"Mz");
    const double* p    = (double*) get_arg(OPT,"P");
    const double* box  = (double*) get_arg(OPT,"box");
    const double* a    = (double*) get_arg(OPT,"a"); /* z-dir offset. RENAME */
#ifndef KAISER
    const double* c    = (double*) get_arg(OPT,"c");
    params->c = *c;
    params->d = pow(params->c/PI,1.5);
#else
    const double* beta = (double*) get_arg(OPT,"beta");
    params->beta = *beta;
    char* window = (char*) get_arg_str(OPT,"window");
    params->window = get_kaiserlike_window(window);
#ifdef FORCE
    params->window_deriv = get_kaiserlike_window_deriv(window);
#endif
    if (strcmp(window, "kaiser_poly") == 0) {
      params->use_polynomial_window = 1;
      params->polynomial_degree = (int) * (double*) get_arg(OPT,"polynomial_degree");
    } else {
      params->use_polynomial_window = 0;
    }
    mxFree(window);
#endif

    params->N = N;
    params->P = (int) *p;
    params->P_half=half( (int) *p );
    params->P_eval_window = params->P;
    params->h = box[0]/m[0];
    params->a = a[0];

    params->dims[0] = (int)  m[0];
    params->dims[1] = (int)  m[0];
    params->dims[2] = (int) mz[0];

    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1]+params->P;
    params->npdims[2] = params->dims[2];

#endif

#ifdef ONE_PERIODIC
    const double* m    = (double*) get_arg(OPT,"M");
    const double* my   = (double*) get_arg(OPT,"My");
    const double* mz   = (double*) get_arg(OPT,"Mz");
    const double* p    = (double*) get_arg(OPT,"P");
    const double* box  = (double*) get_arg(OPT,"box");  
    /* y- and z-dir offsets. */
    const double* a    = (double*) get_arg(OPT,"free_offset");
#ifndef KAISER
    const double* c    = (double*) get_arg(OPT,"c");
    params->c = *c;
    params->d = pow(params->c/PI,1.5);
#else
    const double* beta = (double*) get_arg(OPT,"beta");
    params->beta = *beta;
    char* window = (char*) get_arg_str(OPT,"window");
    params->window = get_kaiserlike_window(window);
#ifdef FORCE
    params->window_deriv = get_kaiserlike_window_deriv(window);
#endif
    if (strcmp(window, "kaiser_poly") == 0) {
      params->use_polynomial_window = 1;
      params->polynomial_degree = (int) * (double*) get_arg(OPT,"polynomial_degree");
    } else {
      params->use_polynomial_window = 0;
    }
    mxFree(window);
#endif

    params->N = N;
    params->P = (int) *p;
    params->P_half=half( (int) *p );
    params->P_eval_window = params->P;
    params->h = box[0]/m[0];
    params->a = a[0];
    params->b = a[1];

    params->dims[0] = (int)  m[0];
    params->dims[1] = (int) my[0];
    params->dims[2] = (int) mz[0];

    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1];
    params->npdims[2] = params->dims[2];

#endif

#ifdef ZERO_PERIODIC
// NB: It seems this is never used; THREE_PERIODIC is used in the 0P case
    const double* m    = (double*) get_arg(OPT,"M");
    const double* p    = (double*) get_arg(OPT,"P");
    const double* box  = (double*) get_arg(OPT,"box");
    const double* a    = (double*) get_arg(OPT,"a"); /* xyz-dir offset.*/
#ifndef KAISER
    const double* c    = (double*) get_arg(OPT,"c");
    params->c = *c;
    params->d = pow(params->c/PI,1.5);
#else
    const double* beta = (double*) get_arg(OPT,"beta");
    params->beta = *beta;
    char* window = (char*) get_arg_str(OPT,"window");
    params->window = get_kaiserlike_window(window);
#ifdef FORCE
    params->window_deriv = get_kaiserlike_window_deriv(window);
#endif
    if (strcmp(window, "kaiser_poly") == 0) {
      params->use_polynomial_window = 1;
      params->polynomial_degree = (int) * (double*) get_arg(OPT,"polynomial_degree");
    } else {
      params->use_polynomial_window = 0;
    }
    mxFree(window);
#endif

    params->N = N;
    params->P = (int) *p;
    params->P_half=half( (int) *p );
    params->P_eval_window = params->P;
    params->h = box[0]/m[0];
    params->a = a[0]; // FIXME: we assume the same offset in each direction!

    params->dims[0] = (int) m[0];
    params->dims[1] = (int) m[1];
    params->dims[2] = (int) m[2];

    params->npdims[0] = params->dims[0];
    params->npdims[1] = params->dims[1];
    params->npdims[2] = params->dims[2];

#endif
}
