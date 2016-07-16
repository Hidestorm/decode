/*
 * File: nxt_stat.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

#ifndef NXT_STAT_H
#define NXT_STAT_H

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
extern void nxt_stat(double current_state, double input, double L, double k,
                     double *next_state, emxArray_real_T *memory_contents);

#endif

/*
 * File trailer for nxt_stat.h
 *
 * [EOF]
 */
