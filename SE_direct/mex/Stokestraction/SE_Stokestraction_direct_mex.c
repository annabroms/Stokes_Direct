#include "mex.h"
#include "SE_Stokestraction_direct.h"
#include "../SE_direct_common.c"

#define X   prhs[0] // Source locations
#define Q   prhs[1] // Source strengths
#define NT  prhs[2] // Target normals
#define OPT prhs[3] // Ewald options

#define U   plhs[0] // Output

#ifndef VERBOSE
#define VERBOSE 0
#endif

// Select what to compute from command-line defined variables.
// Typical compiler defines:
// (...) -DEWALD_REAL -DTHREE_PERIODIC
#ifdef EWALD_FULL
#ifdef ONE_PERIODIC
#define EWALD_KERNEL SE1P_Stokestraction_direct_full
#define EWALD_TAG "1P_FULL"
#elif ZERO_PERIODIC
#define EWALD_KERNEL SE0P_Stokestraction_direct_full
#define EWALD_TAG "0P_FULL"
#else
#error "Must give -D<ZERO/ONE>_PERIODIC to compiler"
#endif
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
    // Check input
    if (nrhs != 4)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "Wrong number of input arguments (x, q, n, opt)");

    // Input dimensions, X and Q should be N-by-3
    if (mxGetNumberOfElements(X) > 0 && mxGetN(X) != 3)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "x must be N-by-3");
    const size_t N = mxGetM(X);
    if (mxGetM(Q) != N || mxGetN(Q) != mxGetN(X))
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "q must be N-by-3");

    // Parameter values
    ewald_opts opt;
    unpack_mex_opt(&opt, OPT);
    // Check options
#ifdef THREE_PERIODIC
    if (opt.box[0] <= 0 || opt.box[1] <= 0 || opt.box[2] <= 0)
#elif TWO_PERIODIC
    if (opt.box[0] <= 0 || opt.box[1] <= 0)
#elif ONE_PERIODIC
    if (opt.box[0] <= 0)
#else
    if (false)
#endif
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.box missing or invalid");
#ifndef EWALD_FULL
    if (opt.xi <= 0)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.xi missing or invalid");
#endif
#ifdef CUTOFF
    if (opt.rc <= 0)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.rc missing or invalid");
#endif
#if !defined(ZERO_PERIODIC) && !defined(EWALD_FOURIER_K0) && !defined(EWALD_SELF)
    if (opt.layers == -1)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "opt.layers missing or invalid");
#endif

    // Pointers
#ifdef EXTERNAL
    size_t num_eval = opt.N_ext;
#else
    size_t num_eval = opt.N_eval;
    if (num_eval == -1) num_eval = N;
#endif

#ifdef EXTERNAL
    if (mxGetM(NT) != num_eval || mxGetN(NT) != 3)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "n must be N_ext-by-3 for external evaluation");
#else
    if (mxGetM(NT) != N || mxGetN(NT) != 3)
        mexErrMsgIdAndTxt("SE:DirectMEXOpt", "n must be N-by-3");
#endif

#ifndef EWALD_SELF
    const double* x = mxGetPr(X);
#endif
    const double* q = mxGetPr(Q);
    const double* n = mxGetPr(NT);
    U = mxCreateDoubleMatrix(num_eval, 3, mxREAL); // output
    double* restrict u = mxGetPr(U);

    if (VERBOSE)
    {
#ifdef EXTERNAL
        __PRINTF("\n[DIRECT STOKESTRACTION EWALD (%s EXT)] MEX N=(%d,%d) ", EWALD_TAG, N, num_eval);
#elif CUTOFF
        __PRINTF("\n[DIRECT STOKESTRACTION EWALD (%s CUT)] MEX N=(%d,%d) ", EWALD_TAG, N, num_eval);
#else
        __PRINTF("\n[DIRECT STOKESTRACTION EWALD (%s)] MEX N=(%d,%d) ", EWALD_TAG, N, num_eval);
#endif
        __PRINTF("xi=%.2f layers=%d\n", opt.xi, opt.layers);
    }

    // Call kernel
#ifdef EWALD_SELF
    EWALD_KERNEL(u, num_eval, q, n, N, &opt);
#else
    EWALD_KERNEL(u, num_eval, x, q, n, N, &opt);
#endif

    free_ewald_opts(&opt);
}
