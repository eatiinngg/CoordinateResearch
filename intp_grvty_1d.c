#include "intp_grvty_1d.h"
#include "main.h"
#include "intp_gauss_1d.h"

const int scale[AXIS] = {116, 82};

int IntpalgGv1D(int * d, const int center, const int axis)
{
    int delta = 0;

    delta = ((d[2] - d[0]) * scale[axis]) / (d[0]+d[2]);

    int tmpPos = delta + (center << INTER_PRECISE);

    return tmpPos;
}

