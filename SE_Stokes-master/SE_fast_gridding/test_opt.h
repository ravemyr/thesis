#include "SE_fgg.h"

void set_params(SE_FGG_params* params, int N, double box[3], int m[3],
	    int p, double c, double beta) {

    params->N = N;
    params->P =  p;
    params->P_half=half(p);
    params->c = c;
    params->d = pow(params->c/PI,1.5);
    params->h = box[0]/m[0];
    params->a = -FGG_INF;
    params->beta = beta;

    params->dims[0] = m[0];
    params->dims[1] = m[1];
    params->dims[2] = m[2];

    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1]+params->P;
    params->npdims[2] = params->dims[2]+params->P;
}
