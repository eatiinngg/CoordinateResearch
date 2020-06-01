#ifndef INTP_GAUSS_1D_H
#define INTP_GAUSS_1D_H

# define GAUSS_PRECISE 15
# define GP (1<<GAUSS_PRECISE)
# define IP (1<<INTER_PRECISE)

#define INTER_PRECISE 7
enum
{
    NLGM_1D_AMP = 0,
    NLGM_1D_MU = 1
};

typedef struct
{
    int max_iter;
    int iter;
    int rs2;                    // 1/sigma^2
    int args[2];          // amp and mu
 //   int pos;          // amp and mu
} nlgm_1d_ctx_t;

#define EXP_SHIFT_X 7
#define EXP_SHIFT_VAL 15

#define EXP_MUL_X 0x80
#define EXP_MUL_VAL 0x8000

#if defined(_USI_DELTA_FROM_TABLE_1D)&&defined(_USI_TILT_PREDICT)
// int IntpalgG1D(int * d, int center, int rs2, int axis, int angle);
int IntpalgG1D(int * d, int center, int rs2, int axis);
#else
int IntpalgG1D(int * d, int center, int rs2);
#endif

#endif /* end of include guard: INTP_GAUSS_1D_H */
