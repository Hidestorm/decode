/*
 * File: vi_de_rtwutil.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

#ifndef VI_DE_RTWUTIL_H
#define VI_DE_RTWUTIL_H

/* Include Files */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "omp.h"
#include "vi_de_types.h"

/* Function Declarations */
extern int div_s32_floor(int numerator, int denominator);
extern double rt_powd_snf(double u0, double u1);
extern double rt_remd_snf(double u0, double u1);
extern double rt_roundd_snf(double u);

#endif

/*
 * File trailer for vi_de_rtwutil.h
 *
 * [EOF]
 */
