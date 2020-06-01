#include <stdio.h>
#include <stdlib.h> ///hi~~
#include <assert.h>
#include <string.h>
#include "main.h"
#include "intp_gauss_1d.h"
#include "intp_grvty_1d.h"
#define SAMPLE_1mm 1
#define ALG_MIX

#define ALG_DIFF_TH_TX 927
#define ALG_DIFF_TH_RX 1122

int flag(int, int, int);

int main(int argc, const char * argv[])
{
    int sigma[AXIS] = {52, 47};
    int Pos[AXIS] = {0, 0};
    int data[9] = {0};

    FILE * fp = fopen("data_tap_a_0304","r");
    /* FILE * fp = fopen("data_tap_a_0304_wired","r"); */
    assert(1);
    char str[STRMAX] = {0};
    char * sub = NULL;

    int frame_cnt = 1;
    int d[3] = {0};
    int index[AXIS] = {0};
    const char * const symbol = ",";
#ifdef ALG_MIX
    int alg_diff_th[AXIS] = {ALG_DIFF_TH_TX, ALG_DIFF_TH_RX};
#endif

#if SAMPLE_1mm
    int j = 0;
    int k = 0;
    int frame_flag = 1;
    int cal_flag = 0;
    int frame_cnt_mm = 0;
#endif
    while ((fgets(str,STRMAX,fp)) != NULL)
    {
        if (str[0] < '0' || str[0] > '9')
        {
            frame_cnt++;
#if SAMPLE_1mm
            frame_flag = 1;
#endif
            continue;
        }
#if SAMPLE_1mm
        if (frame_flag)
        {
            cal_flag = flag(j, k, frame_cnt);
            if (frame_cnt % 2 == 1) j++;
            if (frame_cnt % 201 == 0) j = 0;
            if (frame_cnt % 402 == 0) k++;
            frame_flag = 0;
            frame_cnt_mm += cal_flag;
        }
        if (cal_flag)
        {
            data[0] = frame_cnt_mm;
#endif
            int i = 0;
#if !SAMPLE_1mm
            data[0] = frame_cnt;
            //printf("1\n");
#endif
            sub = strtok(str, symbol);
            while (sub)
            {
                printf("%d", data[i]);
                if (i == 8) break;
                printf(",");
                i++;
                sub = strtok(NULL, symbol);
                data[i] = atoi(sub);
            }
            int num = 0;
            for (int axis = 0; axis < AXIS; axis++)
            {
                num = 4*axis;
                index[axis] = data[num + 1];
                d[0] = data[num + 2];
                d[1] = data[num + 3];
                d[2] = data[num + 4];
#ifdef ALG_MIX
                if (d[1] < alg_diff_th[axis])
                {
                    Pos[axis] = IntpalgG1D(d, index[axis], sigma[axis]);
                }
                else
                {
                    Pos[axis] = IntpalgGv1D(d, index[axis], axis);
                }

#else
                /* Pos[axis] = IntpalgG1D(d, index[axis], sigma[axis]); */
                Pos[axis] = IntpalgGv1D(d, index[axis], axis);
#endif
            }
            printf(",%d,%d\n",Pos[0], Pos[1]);
#if SAMPLE_1mm
        }
#endif
    }
    fclose(fp);

    return 0;
}

int flag(int j,int k,int cnt)
{
    return (402 * k + 2 * j + 1 == cnt);
}
