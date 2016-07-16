/*
 * File: vi_de.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

#ifndef VI_DE_H
#define VI_DE_H

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
extern void vi_de(const double G[6], double k, const creal_T channel_output[1002],
                  const creal_T s16[16], emxArray_real_T *decoder_output);

#endif

/*
 * File trailer for vi_de.h
 *
 * [EOF]
 */
