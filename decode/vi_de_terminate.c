/*
 * File: vi_de_terminate.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "vi_de.h"
#include "vi_de_terminate.h"
#include "vi_de_data.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void vi_de_terminate(void)
{
  omp_destroy_nest_lock(&emlrtNestLockGlobal);
}

/*
 * File trailer for vi_de_terminate.c
 *
 * [EOF]
 */
