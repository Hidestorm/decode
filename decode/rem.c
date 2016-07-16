/*
 * File: rem.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "vi_de.h"
#include "rem.h"
#include "vi_de_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const double x[2]
 *                double r[2]
 * Return Type  : void
 */
void b_rem(const double x[2], double r[2])
{
  int k;
  for (k = 0; k < 2; k++) {
    r[k] = rt_remd_snf(x[k], 2.0);
  }
}

/*
 * File trailer for rem.c
 *
 * [EOF]
 */
