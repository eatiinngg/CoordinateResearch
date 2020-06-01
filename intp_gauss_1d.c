//#include "main.h"
//#include "usi_vote_alg.h"
//#include "usi_slot.h"
//#include "usi_delta_table.h"
#include "intp_gauss_1d.h"

int numerExp(int minus_x)
{
#define TALBE_SIZE 10
// iput:2^7   output:2^15
// goal: approximate exp(x) by Tarlor series
// expand exp(x) at point a

    if (minus_x > 0)
    {
        return -1; //invalid range
    }
    else if (minus_x == 0)
    {
        return 32768;
    }

// exp table: exp(n)*2^15, n = 0,1,...,-9
    static const int etable[TALBE_SIZE] =
    {
        32768,
        12055,
        4435,
        1631,
        600,
        221,
        81,
        30,
        11,
        4
    };

    int h;
    int h_sign;
    int p;

    int x = -minus_x;  // x > 0
// step 0:
// find an integer p such that p-1 < x < p
// denote I = INTER_PRECISE
// the input number minus_x has fixed-point accuracy (INTER_PRECISE)
// therefore we look for bound such that
// bd_left < minus_x < bd_right  (eq 1)
// Hence  p-1 < minus_x/(2^I) < p     (eq 1)

    int p_right = -(x >> EXP_SHIFT_X);
    int p_left = p_right - 1;
    int bd_right = -( -p_right << EXP_SHIFT_X);
    int bd_left  = -( -p_left  << EXP_SHIFT_X);
    if ( -p_right >= TALBE_SIZE)
    {
        //exceed table size!
        return 0;
    }

    //( check bd_r - minus_x >? minus_x - bd_l )
    if (bd_left - 2*minus_x + bd_right > 0)
    {
        h = minus_x - bd_left;
        h_sign = 1;
        p = p_left;
    }
    else
    {
        h = bd_right - minus_x;
        h_sign = -1;
        p = p_right;
    }

    // h = h_sign* h;

    if (h == 0)
        return etable[-p];
    // the Taylor expansion of f(x) is
    // f(x) = f(a) + f'(a)*(x-a) + (f''(a)/2)*(x-a)^2 + (f'''(a)/3!)*(x-a)^3
    // e(x) = e(a) + e(a)*(x-a) + e(a)((x-a)^2)/2 + e(a)((x-a)^3)/3!
    //      = e(a)*(1 + (x-a) + ((x-a)^2)/2 + ((x-a)^3)/3! )

    // dem_exp =`(1 + (x-a) + ((x-a)^2)/2 + ((x-a)^3)/3! )
    int dem_exp
        = EXP_MUL_X
        + h_sign*h
        + (h*h  >> (EXP_SHIFT_X + 1))
        + h_sign*(h*h*h >> (2 *EXP_SHIFT_X + 1 ))/3;

    //      e(a)*(1 + (x-a) + ((x-a)^2)/2 + ((x-a)^3)/3! )
    return (etable[-p] * dem_exp) >> EXP_SHIFT_X;
}

int Nlgm1DCalModle(int *g, int *J, int axis, nlgm_1d_ctx_t * gauss)
{
    int i;
    int u;
    int d;
    int x[3] = { -1, 0, 1};
  //  int y[9] = { -1, -1, -1, 0, 0, 0, 1, 1, 1 };

    // g = A * exp((-(x - xc).^2 - (y - yc).^2)/(2*s^2));
    // jxc = ((x - xc)/(s^2)) .* g;
    // jyc = ((y - yc)/(s^2)) .* g;

    for (i = 0; i < 3; ++i)
    {
        // INTER_PRECISE domain
        d = ((x[i] << INTER_PRECISE) - gauss->args[NLGM_1D_MU]);
     //   d[AXIS_Y] = ((y[i] << INTER_PRECISE) - gauss->pos[AXIS_Y]);

#if 0 //single sigma
        u = (d[0] * d[0] + d[1] * d[1]) >> INTER_PRECISE;
        u = (u * gauss->rs2) >> (GAUSS_PRECISE + 1);
#endif
        u  = ((((d*d)) >> INTER_PRECISE) * gauss->rs2) >> (GAUSS_PRECISE + 1);
   //     u += ((((d[1]*d[1])) >> INTER_PRECISE) * gauss->rs2[1]) >> (GAUSS_PRECISE + 1);


        // INTER_PRECISE -> GAUSS_PRECISE domain
        //u = ( u * gauss->rs2 ) >> ( INTER_PRECISE + 1 );
        //g[i] = expn ( -u, 8 );

        g[i] = numerExp(-u);
        //  printf( "%d ", u );

        if (axis == NLGM_1D_AMP)
        {
            J[i] = g[i];
        }

        g[i] *= gauss->args[NLGM_1D_AMP];
        g[i] >>= GAUSS_PRECISE;

        if (axis == NLGM_1D_MU)
        {
            J[i] = (d * gauss->rs2) >> GAUSS_PRECISE;
            J[i] = (J[i] * g[i]) >> INTER_PRECISE;
        }
    }


    return 0;
}

inline static
int vecDot(int *v1, int *v2)
{
    int i;
    int s = 0;

    for (i = 0; i < 3; ++i)
    {
        s += v1[i] * v2[i];
    }

    return s;
}

int Nlgm1DGaussProc(int *d, nlgm_1d_ctx_t * gauss)
{
    int i, iter;
    int axis;

    int g[3];
    int f[3];
    int J[3];
    int JJ;
    int Jf;
    int dp;
    int dpnorm;
    int tol = 4;

    iter = 0;

    const int delta_range = 63;

    do
    {
        //printf( "\niter %d:\n", iter );
        dpnorm = 0;
        //for ( axis = 0; axis < AXIS_NUM +1; ++axis )
        for (axis = 0; axis < 2; ++axis)
        {
            // calculate Gaussian and Jacobian vectors
            //
            //  g = A * exp((-(x - xc).^2 - (y - yc).^2)/(2*s^2));
            //  J = ((x - xc)/(s^2)) .* g;
            //   or ((y - yc)/(s^2)) .* g;
            //
            Nlgm1DCalModle(g, J, axis, gauss);

            // f = g - d
            for (i = 0; i < 3; ++i)
            {
                f[i] = g[i] - d[i];
            }

            // Gauss-Newton method
            //
            // (J'J)dp = -(J')f
            //
            JJ = vecDot(J, J);
            Jf = vecDot(J, f);

            switch (axis)
            {
                case NLGM_1D_MU:
                  dp = JJ?((-Jf << INTER_PRECISE) / JJ): 0;
                  gauss->args[axis] += dp;

                  int * pos = &gauss->args[NLGM_1D_MU];

                  if (*pos >= delta_range)  *pos  =  delta_range;
                  if (*pos <= -delta_range) *pos = -delta_range;

                  dpnorm += dp * dp;
                  break;

                case NLGM_1D_AMP:
                  dp = JJ?(-Jf / (JJ >> GAUSS_PRECISE)):0;
                  gauss->args[NLGM_1D_AMP] += dp;
                  break;
            }

        }

        iter++;
    }
    while (iter < gauss->max_iter && dpnorm > tol);

    return iter;
}

#if defined(_USI_DELTA_FROM_TABLE_1D)&&defined(_USI_TILT_PREDICT)
/* int IntpalgG1D(int * d, int center, int rs2, int axis, int angle) */
int IntpalgG1D(int * d, int center, int rs2, int axis)
#else
int IntpalgG1D(int * d, int center, int rs2)
#endif
{

    nlgm_1d_ctx_t nlgm_ctx;
    nlgm_1d_ctx_t * pnlgm_ctx = &nlgm_ctx;

    /* GAUSS2D inst_gauss; */
    /* GAUSS2D * gauss = &inst_gauss; */

    /* s is the real std. Usually s is belongs to (0,1].
     * x is the arg for the tuning-setting
     * GP == 2^15, for fixed-point method
     *
     * s = x/100
     * s^2 = x^2/100^2
     * 1/s^2 = 100^2/x^2
     * rs2 = 2^15*100^2/x^2
     */

    //int std[2] = {gx_MutualConfig.Intpalg_nlgm_sx, gx_MutualConfig.Intpalg_nlgm_sy};
    int std = rs2;

    //    intpAlgGauss2DFittingDynamicaStd(std);

    const int ts = std * std;
    //    const int tsy = std[1] * std[1];

    pnlgm_ctx->rs2 = (GP*10000)/ts;

    //gauss->amp = pdif[x0][y0] + 1;
    pnlgm_ctx->args[NLGM_1D_AMP] = d[1] + 1;

    //TODO HSIN
    //setting
    //gauss->max_iter = INTPALG_NLGM_MAX_ITER;
    pnlgm_ctx->max_iter = 20;

    pnlgm_ctx->args[NLGM_1D_MU] = 0;
    pnlgm_ctx->iter = Nlgm1DGaussProc(d, pnlgm_ctx);

    int delta = pnlgm_ctx->args[NLGM_1D_MU];

    //int axis = 0;
    //int idx[2] = {x0, y0};


    /* int tmpPos */
    /*     = delta[axis] + (idx[axis] << INTER_PRECISE ) + gx_TopConfig.BoundaryShift[axis]; */

    //Boundary Shift

#if defined(_USI_DELTA_FROM_TABLE_1D)&&defined(_USI_TILT_PREDICT)
    extern int g_VoteAngle;
    frame_AB_t frame_idx = ( iior( REG_RX_CONFIG2 ) >> 24 == USI_TSISR_PRIOR )
        ? (FRAME_A) : (FRAME_B);
    delta = delta_fun(g_VoteAngle, axis, frame_idx, delta);
#endif

    int tmpPos = delta + (center << INTER_PRECISE );

    return tmpPos;
    /* mutualTP[pctx->idx].coord[1-axis] = tmpPos; */

}
