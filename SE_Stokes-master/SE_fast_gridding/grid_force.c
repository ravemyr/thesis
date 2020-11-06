// -----------------------------------------------------------------------------
void 
SE_FGG_grid_split_SSE_dispatch_force(SE_FGG_work* work, 
				     const SE_state *st,
				     const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 0
    // THIS BYPASSES THE FAST SSE KERNELS.
    // 
    // THEY ARE PLATFORM DEPENDENT, AND THUS MAY NOT WORK OUT OF THE BOX.
    // REMOVE THIS BLOCK ONLY IF YOU ARE FAMILIAR WITH BASIC DEBUGGING, 
    // THE BASICS OF SSE INTRINSICS, AND ARE WILLING TO UNDERSTAND WHERE
    // THE (DATA ALIGNMENT) PRECONDITIONS OF SSE INSTRUCTIONS MAY BREAK
    // IN THE SSE CODE BELOW.
    __DISPATCHER_MSG("[FGG GRID SSE] SSE Disabled\n");
    SE_FGG_grid_split_force(work, st, params);
    return;
#endif

    // if P is odd, or if either increment is odd, fall back on vanilla
    if( is_odd(p) || is_odd(incri) || is_odd(incrj) )
    {
	__DISPATCHER_MSG("[FGG GRID SSE] SSE Abort (PARAMS)\n");
	SE_FGG_grid_split_force(work, st, params);
	return;
    }

#if 0
    // If the work arrays zs or zX are misaligned, fall back on vanilla.
    // These arrays are dynamically allocated, so getting this alignment
    // is really the compilers job! Once you trust it, remove this 
    // check, because the long integer modulus operation is not fast.
    if( ( (unsigned long) work->zs)%16 != 0 || 
	( (unsigned long) work->zx)%16 != 0 || 
	( (unsigned long) work->zy)%16 != 0 ||
	( (unsigned long) work->zz)%16 != 0 )
    {
	__DISPATCHER_MSG("[FGG GRID SSE] SSE Abort (DATA)\n");
	SE_FGG_grid_split_force(work, params);
	return;
    }
#endif
    
    // otherwise the preconditions for SSE codes are satisfied. 
    
    if(p==16)
    {
	// specific for p=16
	__DISPATCHER_MSG("[FGG GRID SSE] P=16\n");
	SE_FGG_grid_split_SSE_P16_force(work, st, params); 
    }
      else if(p%8==0)
    {
	// specific for p divisible by 8
	__DISPATCHER_MSG("[FGG GRID SSE] P unroll 8\n");
	SE_FGG_grid_split_SSE_u8_force(work, st, params); 
    }
    else
    {
	// vanilla SSE code (any even p)
	__DISPATCHER_MSG("[FGG GRID SSE] Vanilla\n");
	SE_FGG_grid_split_SSE_force(work, st, params);
    }
}


// -----------------------------------------------------------------------------
void SE_FGG_grid_split_force(SE_FGG_work* work, 
			     const SE_state* st,
			     const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double cij0,qn;
    int idx0, zidx, idxzz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
    {
	qn = st->q[n];
	idx0 = work->idx[n];

	// inline vanilla loop
	zidx = 0;
	for(i = 0; i<p; i++)
	{
	    for(j = 0; j<p; j++)
	    {
		cij0 = zx[p*n+i]*zy[p*n+j];
		idxzz=p*n;
		for(k = 0; k<p; k++)
		{
		    H[idx0] += zs[zidx]*zz[idxzz]*cij0*qn;
		    idx0++; zidx++; idxzz++;
		}
		idx0 += incrj; 
	    }
	    idx0 += incri; 
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_P16_force(SE_FGG_work* work, 
				     const SE_state* st,
				     const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    double qn;
    int idx, idx_zs, i, j;
    const int incrj = params->npdims[2]-16; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-16);// outer increment

    __m128d rZZ0, rZZ1, rZZ2, rZZ3, rZZ4, rZZ5, rZZ6, rZZ7; 
    __m128d rH0, rH1, rH2, rH3;
    __m128d rC, rZS0;

    for(int n=0; n<N; n++)
    {
	qn = st->q[n];
	idx = work->idx[n];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);

	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

        rZZ0 = _mm_load_pd(zz + n*16     );
        rZZ1 = _mm_load_pd(zz + n*16 + 2 );
        rZZ2 = _mm_load_pd(zz + n*16 + 4 );
        rZZ3 = _mm_load_pd(zz + n*16 + 6 );
        rZZ4 = _mm_load_pd(zz + n*16 + 8 );
        rZZ5 = _mm_load_pd(zz + n*16 + 10);
        rZZ6 = _mm_load_pd(zz + n*16 + 12);
        rZZ7 = _mm_load_pd(zz + n*16 + 14);

	if(idx%2 == 0) // H[idx0] is 16-aligned
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

                    /* 0 - 3 */ 
                    rH0  = _mm_load_pd( H+idx    );
                    rH1  = _mm_load_pd( H+idx + 2);
                    rH2  = _mm_load_pd( H+idx + 4);
                    rH3  = _mm_load_pd( H+idx + 6);

                    rZS0 = _mm_load_pd( zs + idx_zs);
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 2);                   
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 4);                   
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 6);                   
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0));

		    _mm_store_pd(H + idx, rH0);
		    _mm_store_pd(H + idx + 2, rH1);
		    _mm_store_pd(H + idx + 4, rH2);
		    _mm_store_pd(H + idx + 6, rH3);

                    /* 4 - 7*/ 
		    rH0  = _mm_load_pd( H+idx + 8 );
                    rH1  = _mm_load_pd( H+idx + 10);
                    rH2  = _mm_load_pd( H+idx + 12);
                    rH3  = _mm_load_pd( H+idx + 14);

                    rZS0 = _mm_load_pd( zs + idx_zs + 8);
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 10);                   
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 12);                   
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 14);                   
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0));

		    _mm_store_pd(H + idx + 8 , rH0);
		    _mm_store_pd(H + idx + 10, rH1);
		    _mm_store_pd(H + idx + 12, rH2);
		    _mm_store_pd(H + idx + 14, rH3);

		    idx += incrj + 16;
		    idx_zs += 16;
		}
		idx += incri;
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

                    /* 0 - 3 */ 
                    rH0  = _mm_loadu_pd( H+idx    );
                    rH1  = _mm_loadu_pd( H+idx + 2);
                    rH2  = _mm_loadu_pd( H+idx + 4);
                    rH3  = _mm_loadu_pd( H+idx + 6);

		    // if zs does not have 16-byte alignment, this will core.
		    // PLATFORM AND COMPILER DEPENDENT (FIXME)
                    rZS0 = _mm_load_pd( zs + idx_zs);
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 2);                   
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 4);                   
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 6);                   
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0));

		    _mm_storeu_pd(H + idx, rH0);
		    _mm_storeu_pd(H + idx + 2, rH1);
		    _mm_storeu_pd(H + idx + 4, rH2);
		    _mm_storeu_pd(H + idx + 6, rH3);

                    /* 4 - 7*/ 
		    rH0  = _mm_loadu_pd( H+idx + 8 );
                    rH1  = _mm_loadu_pd( H+idx + 10);
                    rH2  = _mm_loadu_pd( H+idx + 12);
                    rH3  = _mm_loadu_pd( H+idx + 14);

                    rZS0 = _mm_load_pd( zs + idx_zs + 8);
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 10);                   
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 12);                   
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 14);                   
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0));

		    _mm_storeu_pd(H + idx + 8 , rH0);
		    _mm_storeu_pd(H + idx + 10, rH1);
		    _mm_storeu_pd(H + idx + 12, rH2);
		    _mm_storeu_pd(H + idx + 14, rH3);

		    idx += incrj + 16;
		    idx_zs += 16;
		}
		idx += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_u8_force(SE_FGG_work* work, 
			            const SE_state* st,
				    const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;
    double qn;
    int idx0, idx_zs, idx_zz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    __m128d rH0, rZZ0, rZS0, rC;
    __m128d rH1, rZZ1, rZS1;
    __m128d rH2, rZZ2, rZS2;
    __m128d rH3, rZZ3, rZS3;

    for(int n=0; n<N; n++)
    {
	qn = st->q[n];
	idx0 = work->idx[n];
	_mm_prefetch( (void*) (H+idx0), _MM_HINT_T0);

	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

	if(idx0%2 == 0) // H[idx0] is 16-aligned
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		    idx_zz=p*n;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm_load_pd( H+idx0     );
			rH1  = _mm_load_pd( H+idx0 + 2 );
			rH2  = _mm_load_pd( H+idx0 + 4 );
			rH3  = _mm_load_pd( H+idx0 + 6 );

			rZZ0 = _mm_load_pd( zz + idx_zz     );
			rZZ1 = _mm_load_pd( zz + idx_zz + 2 );
			rZZ2 = _mm_load_pd( zz + idx_zz + 4 );
			rZZ3 = _mm_load_pd( zz + idx_zz + 6 );

			rZS0 = _mm_load_pd( zs + idx_zs    );
			rZS1 = _mm_load_pd( zs + idx_zs + 2);
			rZS2 = _mm_load_pd( zs + idx_zs + 4);
			rZS3 = _mm_load_pd( zs + idx_zs + 6);

			rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));
			rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1));
			rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2));
			rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3));
			
			_mm_store_pd( H+idx0    , rH0 );
			_mm_store_pd( H+idx0 + 2, rH1 );
			_mm_store_pd( H+idx0 + 4, rH2 );
			_mm_store_pd( H+idx0 + 6, rH3 );

			idx0  +=8;
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx0 += incrj;
		}
		idx0 += incri;
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(i = 0; i<p; i++)
	    {
	    	for(j = 0; j<p; j++)
	    	{
	    	    rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	    	    idx_zz=p*n;
		    
	    	    for(k = 0; k<p; k+=8)
	    	    {
	    		rH0  = _mm_loadu_pd( H+idx0     );
	    		rH1  = _mm_loadu_pd( H+idx0 + 2 );
	    		rH2  = _mm_loadu_pd( H+idx0 + 4 );
	    		rH3  = _mm_loadu_pd( H+idx0 + 6 );

	    		rZZ0 = _mm_load_pd( zz + idx_zz     );
	    		rZZ1 = _mm_load_pd( zz + idx_zz + 2 );
	    		rZZ2 = _mm_load_pd( zz + idx_zz + 4 );
	    		rZZ3 = _mm_load_pd( zz + idx_zz + 6 );

	    		rZS0 = _mm_load_pd( zs + idx_zs    );
	    		rZS1 = _mm_load_pd( zs + idx_zs + 2);
	    		rZS2 = _mm_load_pd( zs + idx_zs + 4);
	    		rZS3 = _mm_load_pd( zs + idx_zs + 6);

			rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));
			rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1));
			rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2));
			rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3));

	    		_mm_storeu_pd( H+idx0    , rH0 );
	    		_mm_storeu_pd( H+idx0 + 2, rH1 );
	    		_mm_storeu_pd( H+idx0 + 4, rH2 );
	    		_mm_storeu_pd( H+idx0 + 6, rH3 );

	    		idx0  +=8;
	    		idx_zs+=8;
	    		idx_zz+=8;
	    	    }
	    	    idx0 += incrj;
	    	}
	    	idx0 += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_force(SE_FGG_work* work, 
			     	 const SE_state* st,
				 const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;
    double qn;
    int idx0, idx_zs, idx_zz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    __m128d rH0, rZZ0, rZS0, rC;

    for(int n=0; n<N; n++)
    {
	qn = st->q[n];
	idx0 = work->idx[n];
	idx_zs = 0;

	if(idx0%2 == 0) // H[idx0] is 16-aligned
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		    idx_zz=p*n;
		    for(k = 0; k<p; k+=2)
		    {
			rH0  = _mm_load_pd( H+idx0     );
			rZZ0 = _mm_load_pd( zz + idx_zz     );
			rZS0 = _mm_load_pd( zs + idx_zs    );

			rZZ0 = _mm_mul_pd(rZZ0,rC);
			rZZ0 = _mm_mul_pd(rZZ0,rZS0);
			rH0  = _mm_add_pd(rH0,rZZ0);

			_mm_store_pd( H+idx0    , rH0 );

			idx0  +=2;
			idx_zs+=2; 
			idx_zz+=2;
		    }
		    idx0 += incrj; 
		}
		idx0 += incri; 
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(i = 0; i<p; i++)
	    {
	    	for(j = 0; j<p; j++)
	    	{
	    	    rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	    	    idx_zz=p*n;
	    	    for(k = 0; k<p; k+=2)
	    	    {
	    		rH0  = _mm_loadu_pd( H+idx0 );
	    		rZZ0 = _mm_load_pd( zz + idx_zz );
	    		rZS0 = _mm_load_pd( zs + idx_zs );
	    		rZZ0 = _mm_mul_pd(rZZ0,rC);
	    		rZZ0 = _mm_mul_pd(rZZ0,rZS0);
	    		rH0  = _mm_add_pd(rH0,rZZ0);
	    		_mm_storeu_pd( H+idx0, rH0 );

	    		idx0  +=2;
	    		idx_zs+=2;
	    		idx_zz+=2;
	    	    }
	    	    idx0 += incrj;
	    	}
	    	idx0 += incri;
	    }
	}
    }
}


// -----------------------------------------------------------------------------
#ifdef __AVX__
void 
SE_FGG_grid_split_AVX_dispatch_force(SE_FGG_work* work, 
				     const SE_state *st,
				     const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 0
    // THIS BYPASSES THE FAST AVX KERNELS.
    __DISPATCHER_MSG("[FGG GRID AVX] AVX Disabled\n");
    SE_FGG_grid_split_force(work, st, params);
    return;
#endif

    // if P, or either increments are not divisible by 4, fall back on vanilla
    if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) )
    {
	__DISPATCHER_MSG("[FGG GRID AVX] AVX Abort (PARAMS)\n");
	SE_FGG_grid_split_force(work, st, params);
	return;
    }

#if 0
    // If the work arrays zs or zx are misaligned, fall back on vanilla.
    // These arrays are dynamically allocated, so getting this alignment
    // is really the compilers job! Once you trust it, remove this 
    // check, because the long integer modulus operation is not fast.
    if( ( (unsigned long) work->zs)%32 != 0 || 
	( (unsigned long) work->zx)%32 != 0 || 
	( (unsigned long) work->zy)%32 != 0 ||
	( (unsigned long) work->zz)%32 != 0 )
    {
	__DISPATCHER_MSG("[FGG GRID AVX] AVX Abort (DATA)\n");
	SE_FGG_grid_split_force(work, params);
	return;
    }
#endif
    
    // otherwise the preconditions for AVX codes are satisfied. 
    
    if(p==16)
    {
	// specific for p=16
	__DISPATCHER_MSG("[FGG GRID AVX] P=16\n");
	SE_FGG_grid_split_AVX_P16_force(work, st, params); 
    }
    else if(p==8)
    {
	// specific for p=8
	__DISPATCHER_MSG("[FGG GRID AVX] P=8\n");
	SE_FGG_grid_split_AVX_P8_force(work, st, params); 
    }
    else if(p%8==0)
    {
	// specific for p divisible by 8
	__DISPATCHER_MSG("[FGG GRID AVX] P unroll 8\n");
	SE_FGG_grid_split_AVX_u8_force(work, st, params); 
    }
    else if(p%4==0)
    {
      // specific for p divisible by 4
	__DISPATCHER_MSG("[FGG GRID AVX] P unroll 4\n");
	SE_FGG_grid_split_AVX_force(work, st, params); 
    }
    else
    {
	// vanilla SSE code (any even p)
	__DISPATCHER_MSG("[FGG GRID AVX] Vanilla\n");
	SE_FGG_grid_split_SSE_force(work, st, params);
    }
}


// -----------------------------------------------------------------------------
void SE_FGG_grid_split_AVX_P16_force(SE_FGG_work* work, 
				     const SE_state* st,
				     const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;


    const int incrj = params->npdims[2]-16; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-16);// outer increment

    double qn;
    int idx, idx_zs, i, j;

    __m256d rZZ0, rZZ1, rZZ2, rZZ3; 
    __m256d rH0, rH1, rH2, rH3;
    __m256d rC, rZS0,rZS1,rZS2,rZS3;

    for(int n=0; n<N; n++)
    {
        qn = st->q[n];
	idx = work->idx[n];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

        rZZ0 = _mm256_load_pd(zz + n*16     );
        rZZ1 = _mm256_load_pd(zz + n*16 + 4 );
        rZZ2 = _mm256_load_pd(zz + n*16 + 8 );
        rZZ3 = _mm256_load_pd(zz + n*16 + 12);

	if(idx%4 == 0) // H[idx0] is 32-aligned
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

                    rH0  = _mm256_load_pd( H+idx    );
                    rH1  = _mm256_load_pd( H+idx + 4);
                    rH2  = _mm256_load_pd( H+idx + 8);
                    rH3  = _mm256_load_pd( H+idx + 12);

                    rZS0 = _mm256_load_pd( zs + idx_zs);
		    rZS1 = _mm256_load_pd( zs + idx_zs + 4);
                    rZS2 = _mm256_load_pd( zs + idx_zs + 8);   
                    rZS3 = _mm256_load_pd( zs + idx_zs + 12);    

		    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
		    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		    rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2));
		    rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3));

		    _mm256_store_pd(H + idx,      rH0);
		    _mm256_store_pd(H + idx + 4,  rH1);
		    _mm256_store_pd(H + idx + 8,  rH2);
		    _mm256_store_pd(H + idx + 12, rH3);

		    idx += incrj + 16;
		    idx_zs += 16;
		}
		idx += incri;
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

                    rH0  = _mm256_loadu_pd( H+idx     );
                    rH1  = _mm256_loadu_pd( H+idx + 4 );
                    rH2  = _mm256_loadu_pd( H+idx + 8 );
                    rH3  = _mm256_loadu_pd( H+idx + 12);

                    rZS0 = _mm256_load_pd( zs + idx_zs     );
                    rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
                    rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
                    rZS3 = _mm256_load_pd( zs + idx_zs + 12);

		    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
		    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		    rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2));
		    rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3));

		    _mm256_storeu_pd(H + idx,      rH0);
		    _mm256_storeu_pd(H + idx + 4,  rH1);
		    _mm256_storeu_pd(H + idx + 8,  rH2);
		    _mm256_storeu_pd(H + idx + 12, rH3);

		    idx += incrj + 16;
		    idx_zs += 16;
		}
		idx += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_AVX_P8_force(SE_FGG_work* work, 
				     const SE_state* st,
				     const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    double qn;
    int idx, idx_zs, i, j;
    const int incrj = params->npdims[2]-8; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-8);// outer increment

    __m256d rZZ0, rZZ1; 
    __m256d rH0, rH1;
    __m256d rC, rZS0,rZS1;

    for(int n=0; n<N; n++)
    {
	qn = st->q[n];
	idx = work->idx[n];
	idx_zs = 0;

        rZZ0 = _mm256_load_pd(zz + n*8     );
        rZZ1 = _mm256_load_pd(zz + n*8 + 4 );

	if(idx%4 == 0) // H[idx0] is 32-aligned
	{
	    for(i = 0; i<8; i++)
	    {
		for(j = 0; j<8; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[8*n+i]*zy[8*n+j] );

                    rH0  = _mm256_load_pd( H+idx    );
                    rH1  = _mm256_load_pd( H+idx + 4);

                    rZS0 = _mm256_load_pd( zs + idx_zs);
		    rZS1 = _mm256_load_pd( zs + idx_zs + 4);

		    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
		    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));

		    _mm256_store_pd(H + idx,      rH0);
		    _mm256_store_pd(H + idx + 4,  rH1);

		    idx += incrj + 8;
		    idx_zs += 8;
		}
		idx += incri;
	    }
	}
	else // H[idx0] is 16-aligned, preventing nice vectorization
	{
	    for(i = 0; i<8; i++)
	    {
		for(j = 0; j<8; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[8*n+i]*zy[8*n+j] );

                    rH0  = _mm256_loadu_pd( H+idx     );
                    rH1  = _mm256_loadu_pd( H+idx + 4 );

                    rZS0 = _mm256_load_pd( zs + idx_zs     );
                    rZS1 = _mm256_load_pd( zs + idx_zs + 4 );

		    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
		    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));

		    _mm256_storeu_pd(H + idx,      rH0);
		    _mm256_storeu_pd(H + idx + 4,  rH1);

		    idx += incrj + 8;
		    idx_zs += 8;
		}
		idx += incri;
	    }
	}
    }
}


// -----------------------------------------------------------------------------
void SE_FGG_grid_split_AVX_u8_force(SE_FGG_work* work, 
			            const SE_state* st,
				    const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    double qn;
    int idx0, idx_zs, idx_zz, i, j, k;

    __m256d rH0, rZZ0, rZS0, rC;
    __m256d rH1, rZZ1, rZS1;

    for(int n=0; n<N; n++)
    {
	qn = st->q[n];
	idx0 = work->idx[n];
	_mm_prefetch( (void*) (H+idx0), _MM_HINT_T0); 
	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

	if(idx0%4 == 0) // H[idx0] is 32-aligned
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		    idx_zz=p*n;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm256_load_pd( H+idx0     );
			rH1  = _mm256_load_pd( H+idx0 + 4 );

			rZZ0 = _mm256_load_pd( zz + idx_zz     );
			rZZ1 = _mm256_load_pd( zz + idx_zz + 4 );

			rZS0 = _mm256_load_pd( zs + idx_zs    );
			rZS1 = _mm256_load_pd( zs + idx_zs + 4);

			rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
			rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
			
			_mm256_store_pd( H+idx0    , rH0 );
			_mm256_store_pd( H+idx0 + 4, rH1 );

			idx0  +=8;
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx0 += incrj;
		}
		idx0 += incri;
	    }
	}
	else // H[idx0] is 16-aligned, preventing nice vectorization
	{
	    for(i = 0; i<p; i++)
	    {
	    	for(j = 0; j<p; j++)
	    	{
	    	    rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	    	    idx_zz=p*n;
		    
	    	    for(k = 0; k<p; k+=8)
	    	    {
	    		rH0  = _mm256_loadu_pd( H+idx0     );
	    		rH1  = _mm256_loadu_pd( H+idx0 + 4 );

	    		rZZ0 = _mm256_load_pd( zz + idx_zz     );
	    		rZZ1 = _mm256_load_pd( zz + idx_zz + 4 );

	    		rZS0 = _mm256_load_pd( zs + idx_zs    );
	    		rZS1 = _mm256_load_pd( zs + idx_zs + 4);

			rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
			rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));

	    		_mm256_storeu_pd( H+idx0    , rH0 );
	    		_mm256_storeu_pd( H+idx0 + 4, rH1 );

	    		idx0  +=8;
	    		idx_zs+=8;
	    		idx_zz+=8;
	    	    }
	    	    idx0 += incrj;
	    	}
	    	idx0 += incri;
	    }
	}
    }
}


// -----------------------------------------------------------------------------
void SE_FGG_grid_split_AVX_force(SE_FGG_work* work, 
			     	 const SE_state* st,
				 const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;
    double qn;
    int idx0, idx_zs, idx_zz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    __m256d rH0, rZZ0, rZS0, rC;

    for(int n=0; n<N; n++)
    {
	qn = st->q[n];
	idx0 = work->idx[n];
	idx_zs = 0;

	if(idx0%4 == 0) // H[idx0] is 32-aligned
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		    idx_zz=p*n;
		    for(k = 0; k<p; k+=4)
		    {
			rH0  = _mm256_load_pd( H+idx0     );
			rZZ0 = _mm256_load_pd( zz + idx_zz     );
			rZS0 = _mm256_load_pd( zs + idx_zs    );

			rZZ0 = _mm256_mul_pd(rZZ0,rC);
			rZZ0 = _mm256_mul_pd(rZZ0,rZS0);
			rH0  = _mm256_add_pd(rH0,rZZ0);

			_mm256_store_pd( H+idx0    , rH0 );

			idx0  +=4;
			idx_zs+=4; 
			idx_zz+=4;
		    }
		    idx0 += incrj; 
		}
		idx0 += incri; 
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(i = 0; i<p; i++)
	    {
	    	for(j = 0; j<p; j++)
	    	{
	    	    rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	    	    idx_zz=p*n;
	    	    for(k = 0; k<p; k+=4)
	    	    {
	    		rH0  = _mm256_loadu_pd( H+idx0 );

	    		rZZ0 = _mm256_load_pd( zz + idx_zz );
	    		rZS0 = _mm256_load_pd( zs + idx_zs );

	    		rZZ0 = _mm256_mul_pd(rZZ0,rC);
	    		rZZ0 = _mm256_mul_pd(rZZ0,rZS0);

	    		rH0  = _mm256_add_pd(rH0,rZZ0);
	    		_mm256_storeu_pd( H+idx0, rH0 );

	    		idx0  +=4;
	    		idx_zs+=4;
	    		idx_zz+=4;
	    	    }
	    	    idx0 += incrj;
	    	}
	    	idx0 += incri;
	    }
	}
    }
}
#endif //AVX
