/*
 * File: deci2bin.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "vi_de.h"
#include "deci2bin.h"
#include "rem.h"
#include "vi_de_emxutil.h"
#include "vi_de_rtwutil.h"

/* Function Definitions */

/*
 * 十进制数x转化为二进制数，二进制数至少表示为t位
 * Arguments    : double x
 *                double y[2]
 * Return Type  : void
 */
void b_deci2bin(double x, double y[2])
{
  int i5;
  int i;
  double b_y[2];
  for (i5 = 0; i5 < 2; i5++) {
    y[i5] = 0.0;
  }

  i = 0;
  while ((x >= 0.0) && (i + 1 <= 2)) {
    y[i] = rt_remd_snf(x, 2.0);
    x = (x - y[i]) / 2.0;
    i++;
  }

  for (i5 = 0; i5 < 2; i5++) {
    b_y[i5] = y[1 - i5];
  }

  for (i5 = 0; i5 < 2; i5++) {
    y[i5] = b_y[i5];
  }
}

/*
 * 十进制数x转化为二进制数，二进制数至少表示为t位
 * Arguments    : double x
 *                double t
 *                emxArray_real_T *y
 * Return Type  : void
 */
void deci2bin(double x, double t, emxArray_real_T *y)
{
  int i2;
  int loop_ub;
  unsigned int i;
  int i3;
  int i4;
  emxArray_real_T *b_y;
  i2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)t;
  emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
  loop_ub = (int)t;
  for (i2 = 0; i2 < loop_ub; i2++) {
    y->data[i2] = 0.0;
  }

  i = 1U;
  while ((x >= 0.0) && (i <= t)) {
    y->data[(int)i - 1] = rt_remd_snf(x, 2.0);
    x = (x - y->data[(int)i - 1]) / 2.0;
    i++;
  }

  if (1.0 > t) {
    i2 = 1;
    i3 = 1;
    i4 = 0;
  } else {
    i2 = (int)t;
    i3 = -1;
    i4 = 1;
  }

  emxInit_real_T(&b_y, 2);
  loop_ub = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = div_s32_floor(i4 - i2, i3) + 1;
  emxEnsureCapacity((emxArray__common *)b_y, loop_ub, (int)sizeof(double));
  loop_ub = div_s32_floor(i4 - i2, i3);
  for (i4 = 0; i4 <= loop_ub; i4++) {
    b_y->data[b_y->size[0] * i4] = y->data[(i2 + i3 * i4) - 1];
  }

  i2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
  loop_ub = b_y->size[1];
  for (i2 = 0; i2 < loop_ub; i2++) {
    y->data[y->size[0] * i2] = b_y->data[b_y->size[0] * i2];
  }

  emxFree_real_T(&b_y);
}

/*
 * File trailer for deci2bin.c
 *
 * [EOF]
 */
