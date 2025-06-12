#include <math.h>
#include "SE_Stokeslet_direct.h"


inline double dot(double * a, double * b)
{
    return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

#ifdef EWALD_FULL
/* Equivalent MATLAB code (not external):
M = numel(opt.eval_idx);
N = size(x, 1);
u = zeros(M, 3);
for i=1:M
  target = opt.eval_idx(i);
  source = [1:target-1 target+1:N];
  rvec = bsxfun(@minus, x(target,:), x(source,:));
  ri = 1./sqrt(sum(rvec.^2, 2));
  rdotf = sum(rvec .* f(source,:), 2);
  u(i,:) = sum( ...
               bsxfun(@times, f(source,:), ri) ...
               + bsxfun(@times, rvec, rdotf.*ri.^3), ...
               1);
end
 */
void SE0P_Stokeslet_direct_full(double* restrict u, const size_t nout,
                                const double* restrict x, const double* restrict f,
                                const size_t N, const ewald_opts* restrict opt)
{
#ifdef EXTERNAL
    const double* restrict xt = opt->eval_ext_x;
#else
    const size_t* restrict idx = opt->eval_idx;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t mi=0; mi<nout; mi++)                // for all evaluation points
    {
#ifndef EXTERNAL
        size_t m;
        if (opt->N_eval == -1)
            m = mi;               // default is to evaluate at all source points
        else
            m = idx[mi];                   // indirect indexing OK in outer loop
#endif
        double p[] = {0,0,0};
        for (size_t n=0; n<N; n++)                          // for all particles
        {
#ifndef EXTERNAL
            if (m == n) continue;                                   // skip self
#endif

            double fn[] = {f[n], f[n+N], f[n+2*N]};
#ifdef EXTERNAL
            double  r[] = {xt[mi]-x[n], xt[mi+nout]-x[n+N], xt[mi+2*nout]-x[n+2*N]};
#else
            double  r[] = {x[m]-x[n], x[m+N]-x[n+N], x[m+2*N]-x[n+2*N]};
#endif
            double rdotf = dot(r,fn);
            double ri = 1.0/sqrt(dot(r,r));
            double ri3 = ri*ri*ri;
            p[0] += ri*fn[0] + rdotf*ri3*r[0];
            p[1] += ri*fn[1] + rdotf*ri3*r[1];
            p[2] += ri*fn[2] + rdotf*ri3*r[2];
        }
        u[mi       ] = p[0];
        u[mi+nout  ] = p[1];
        u[mi+2*nout] = p[2];
    }
}
#endif // EWALD_FULL


