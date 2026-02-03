#include <math.h>
#include "SE_Stokestraction_direct.h"

static inline double dot3(const double* a, const double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
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
  ri5 = 1./sqrt(sum(rvec.^2, 2)).^5;
  rdotf = sum(rvec .* f(source,:), 2);
  rdotn = sum(rvec .* n(target,:), 2);
  u(i,:) =  6 * sum( bsxfun(@times, rvec, rdotn.*rdotf.*ri5), 1);
end
 */
void SE0P_Stokestraction_direct_full(double* restrict u, const size_t nout,
                                     const double* restrict x,
                                     const double* restrict f, const double* restrict nt,
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
    for (size_t m=0; m<nout; m++)                   // for all evaluation points
    {
        double um[] = {0, 0, 0};
#ifdef EXTERNAL
        double xm[] = {xt[m], xt[m+nout], xt[m+2*nout]};
        double nm[] = {nt[m], nt[m+nout], nt[m+2*nout]};
#else
        size_t m_idx;
        if (opt->N_eval == -1)
            m_idx = m;            // default is to evaluate at all source points
        else
            m_idx = idx[m];                // indirect indexing OK in outer loop
        double xm[] = {x[m_idx], x[m_idx+N], x[m_idx+2*N]};
        double nm[] = {nt[m_idx], nt[m_idx+N], nt[m_idx+2*N]};
#endif

        double r[3], fn[3];
        for (size_t n=0; n<N; n++)                          // for all particles
        {
#ifndef EXTERNAL
            if (n == m_idx) continue;                               // skip self
#endif
            // EXTERNAL: Assuming that r != 0 in home box

            r[0] = xm[0] - x[n    ];
            r[1] = xm[1] - x[n+N  ];
            r[2] = xm[2] - x[n+2*N];
            fn[0] = f[n    ];
            fn[1] = f[n+N  ];
            fn[2] = f[n+2*N];

            double rdotf = dot3(r, fn);
            double rdotn = dot3(r, nm);
            double ri2 = 1.0 / dot3(r, r);
            double ri = sqrt(ri2);
            double ri5 = ri2*ri2*ri;
            double C = 6.0 * rdotf * rdotn * ri5;

            um[0] += C*r[0];
            um[1] += C*r[1];
            um[2] += C*r[2];
        }
        u[m       ] = um[0];
        u[m+nout  ] = um[1];
        u[m+2*nout] = um[2];
    }
}
#endif // EWALD_FULL
