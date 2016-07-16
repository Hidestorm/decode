/*
 * File: nxt_stat.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "vi_de.h"
#include "nxt_stat.h"
#include "vi_de_emxutil.h"
#include "deci2bin.h"
#include "vi_de_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : double current_state
 *                double input
 *                double L
 *                double k
 *                double *next_state
 *                emxArray_real_T *memory_contents
 * Return Type  : void
 */
void nxt_stat(double current_state, double input, double L, double k, double
              *next_state, emxArray_real_T *memory_contents)
{
  emxArray_real_T *binary_state;
  emxArray_real_T *binary_input;
  double bnew;
  int apnd;
  emxArray_real_T *next_state_binary;
  int i1;
  int absa;
  int anew;
  int ndbl;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  unsigned int uv0[2];
  int b_k;
  int c_k;
  emxArray_real_T *b;
  emxInit_real_T(&binary_state, 2);
  emxInit_real_T(&binary_input, 2);
  deci2bin(current_state, k * (L - 1.0), binary_state);
  deci2bin(input, k, binary_input);
  bnew = (L - 2.0) * k;
  if (1.0 > bnew) {
    apnd = 0;
  } else {
    apnd = (int)bnew;
  }

  emxInit_real_T(&next_state_binary, 2);
  i1 = next_state_binary->size[0] * next_state_binary->size[1];
  next_state_binary->size[0] = 1;
  next_state_binary->size[1] = binary_input->size[1] + apnd;
  emxEnsureCapacity((emxArray__common *)next_state_binary, i1, (int)sizeof
                    (double));
  absa = binary_input->size[1];
  for (i1 = 0; i1 < absa; i1++) {
    next_state_binary->data[next_state_binary->size[0] * i1] =
      binary_input->data[binary_input->size[0] * i1];
  }

  for (i1 = 0; i1 < apnd; i1++) {
    next_state_binary->data[next_state_binary->size[0] * (i1 +
      binary_input->size[1])] = binary_state->data[i1];
  }

  /*  将二进制数转化为十进制数 */
  if (next_state_binary->size[1] - 1 < 0) {
    ndbl = 0;
    anew = -1;
    bnew = 0.0;
  } else {
    anew = next_state_binary->size[1] - 1;
    ndbl = (int)floor((0.0 - ((double)next_state_binary->size[1] - 1.0)) / -1.0
                      + 0.5);
    apnd = (next_state_binary->size[1] - ndbl) - 1;
    absa = (int)fabs((double)next_state_binary->size[1] - 1.0);
    if (absa >= 0) {
    } else {
      absa = 0;
    }

    if (fabs(0.0 - (double)apnd) < 4.4408920985006262E-16 * (double)absa) {
      ndbl++;
      bnew = 0.0;
    } else if (-apnd > 0) {
      bnew = ((double)next_state_binary->size[1] - 1.0) + -((double)ndbl - 1.0);
    } else {
      ndbl++;
      bnew = apnd;
    }
  }

  emxInit_real_T(&y, 2);
  i1 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)y, i1, (int)sizeof(double));
  if (ndbl > 0) {
    y->data[0] = anew;
    if (ndbl > 1) {
      y->data[ndbl - 1] = bnew;
      apnd = (ndbl - 1) / 2;
      for (absa = 1; absa < apnd; absa++) {
        y->data[absa] = anew - absa;
        y->data[(ndbl - absa) - 1] = bnew - (-(double)absa);
      }

      if (apnd << 1 == ndbl - 1) {
        y->data[apnd] = ((double)anew + bnew) / 2.0;
      } else {
        y->data[apnd] = anew - apnd;
        y->data[apnd + 1] = bnew - (-(double)apnd);
      }
    }
  }

  emxInit_real_T(&b_y, 2);
  i1 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)b_y, i1, (int)sizeof(double));
  apnd = y->size[0] * y->size[1];
  for (i1 = 0; i1 < apnd; i1++) {
    b_y->data[i1] = y->data[i1];
  }

  for (i1 = 0; i1 < 2; i1++) {
    uv0[i1] = (unsigned int)y->size[i1];
  }

  i1 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)y, i1, (int)sizeof(double));
  apnd = b_y->size[1];

#pragma omp parallel for \
 num_threads(omp_get_max_threads()) \
 private(c_k)

  for (b_k = 1; b_k <= apnd; b_k++) {
    c_k = b_k;
    y->data[c_k - 1] = rt_powd_snf(2.0, b_y->data[c_k - 1]);
  }

  emxFree_real_T(&b_y);
  emxInit_real_T1(&b, 1);
  i1 = b->size[0];
  b->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)b, i1, (int)sizeof(double));
  apnd = y->size[1];
  for (i1 = 0; i1 < apnd; i1++) {
    b->data[i1] = y->data[y->size[0] * i1];
  }

  emxFree_real_T(&y);
  if ((next_state_binary->size[1] == 1) || (b->size[0] == 1)) {
    bnew = 0.0;
    for (i1 = 0; i1 < next_state_binary->size[1]; i1++) {
      bnew += next_state_binary->data[next_state_binary->size[0] * i1] * b->
        data[i1];
    }

    *next_state = bnew;
  } else {
    bnew = 0.0;
    for (i1 = 0; i1 < next_state_binary->size[1]; i1++) {
      bnew += next_state_binary->data[next_state_binary->size[0] * i1] * b->
        data[i1];
    }

    *next_state = bnew;
  }

  emxFree_real_T(&b);
  emxFree_real_T(&next_state_binary);
  i1 = memory_contents->size[0] * memory_contents->size[1];
  memory_contents->size[0] = 1;
  memory_contents->size[1] = binary_input->size[1] + binary_state->size[1];
  emxEnsureCapacity((emxArray__common *)memory_contents, i1, (int)sizeof(double));
  apnd = binary_input->size[1];
  for (i1 = 0; i1 < apnd; i1++) {
    memory_contents->data[memory_contents->size[0] * i1] = binary_input->
      data[binary_input->size[0] * i1];
  }

  apnd = binary_state->size[1];
  for (i1 = 0; i1 < apnd; i1++) {
    memory_contents->data[memory_contents->size[0] * (i1 + binary_input->size[1])]
      = binary_state->data[binary_state->size[0] * i1];
  }

  emxFree_real_T(&binary_input);
  emxFree_real_T(&binary_state);
}

/*
 * File trailer for nxt_stat.c
 *
 * [EOF]
 */
