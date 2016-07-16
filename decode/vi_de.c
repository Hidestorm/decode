/*
 * File: vi_de.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "vi_de.h"
#include "rem.h"
#include "nxt_stat.h"
#include "deci2bin.h"
#include "vi_de_emxutil.h"
#include "vi_de_rtwutil.h"

/* Function Definitions */

/*
 * 维特比译码实现程序
 * Arguments    : const double G[6]
 *                double k
 *                const creal_T channel_output[1002]
 *                const creal_T s16[16]
 *                emxArray_real_T *decoder_output
 * Return Type  : void
 */
void vi_de(const double G[6], double k, const creal_T channel_output[1002],
           const creal_T s16[16], emxArray_real_T *decoder_output)
{
  emxArray_real_T *input;
  double L;
  double number_of_states;
  int i0;
  int ic;
  emxArray_real_T *nextstate;
  emxArray_real_T *output;
  int n;
  emxArray_real_T *dec_output_bin;
  emxArray_real_T *state_metric;
  double d0;
  int l;
  int b_l;
  double next_state;
  emxArray_real_T *survivor_state;
  double b[6];
  int b_k;
  double y[2];
  int nx;
  double mtmp;
  int br;
  double binary_output[2];
  emxArray_real_T *survivor_first_bit;
  int ar;
  int ib;
  int idx;
  emxArray_real_T *survivor_second_bit;
  int m;
  emxArray_real_T *b_state_metric;
  emxArray_real_T *c_state_metric;
  double b_m;
  double step;
  double state_sequence[1003];
  double b_n;
  double first_bit_sequence[1002];
  double second_bit_sequence[1002];
  double parallel_binary_output[4];
  emxArray_real_T *decoder_output_matrix;
  double parallel_binary_output_two[4];
  double parallel_binary_output_three[4];
  static const signed char b_b[4] = { 8, 4, 2, 1 };

  double parallel_binary_output_four[4];
  double a;
  double b_a;
  double branch_metric;
  double c_a;
  double d_a;
  boolean_T empty_non_axis_sizes;
  double d1;
  double first_bit_sequence_data[1002];
  double second_bit_sequence_data[1002];
  double e_a;
  double f_a;
  double g_a;
  double parallel_diatance_output[4];
  boolean_T exitg2;
  signed char ii_data[4];
  int ii_size[2];
  static const signed char iv0[2] = { 1, 4 };

  boolean_T exitg1;
  boolean_T guard1 = false;
  double b_parallel_binary_output[16];
  double parallel_double_bit_output_data[16];
  emxInit_real_T(&input, 2);

  /* VITERBI _decorder   卷积码的维特比软判决解码器  */
  /* 输入参数：G：生成矩阵，G是一个n×LK矩阵，该矩阵的每一行确定了从移位计错器到第n个输出间的连接 */
  /*           k：并行输入比特数 */
  /*           channel_output：需要解码的数据 */
  /* 输出参数：decoder_output：解码数据输出（已经截掉尾部的零） */
  /* survivor_state为记录最佳路径的矩阵,cumulated_metric为累积最小汉明距 */
  /* 程序中调用的外部函数：nxt_stat.m：利用当前状态和输入数据求网格图的下一转移状态 */
  /*                       bin2deci.m：二进制转十进制 */
  /*                       deci2bin.m：十进制转二进制 */
  /*                       metric.m  ：计算欧式距离 */
  /* 该译码思想是按照最小欧式距准则确定幸存路径，然后根据幸存路径由后到前得到各级的状态，最后根据各级状态和状态转移图来确定各级的输入 */
  /* 算出n的大小 */
  /*  检查输入各矩阵格式是否合格 */
  /* 确定L，用十进制数表示寄存器状态 */
  L = 3.0 / k;
  number_of_states = rt_powd_snf(2.0, (L - 1.0) * k);

  /* 产生网格图 */
  /* 产生状态转移矩阵，输出矩阵和输入矩阵 */
  /* 在nextstate(状态数x2)(即状态转移矩阵)中， 两列分别为i状态下输入0，1所得到的下一状态（为十进制数） */
  /* 在input(状态数x状态数)矩阵中保存状态的转移关系，input(i,j)是状态i到状态j的输入 */
  /* 在output(状态数x2)(即输出矩阵)中，两列分别为i状态下输入0，1理想输出（为十进制数） */
  i0 = input->size[0] * input->size[1];
  input->size[0] = (int)number_of_states;
  input->size[1] = (int)number_of_states;
  emxEnsureCapacity((emxArray__common *)input, i0, (int)sizeof(double));
  ic = (int)number_of_states * (int)number_of_states;
  for (i0 = 0; i0 < ic; i0++) {
    input->data[i0] = 0.0;
  }

  emxInit_real_T(&nextstate, 2);
  i0 = nextstate->size[0] * nextstate->size[1];
  nextstate->size[0] = (int)number_of_states;
  nextstate->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)nextstate, i0, (int)sizeof(double));
  ic = (int)number_of_states << 1;
  for (i0 = 0; i0 < ic; i0++) {
    nextstate->data[i0] = 0.0;
  }

  emxInit_real_T(&output, 2);
  i0 = output->size[0] * output->size[1];
  output->size[0] = (int)number_of_states;
  output->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)output, i0, (int)sizeof(double));
  ic = (int)number_of_states << 1;
  for (i0 = 0; i0 < ic; i0++) {
    output->data[i0] = 0.0;
  }

  n = 0;
  emxInit_real_T(&dec_output_bin, 2);
  while (n <= (int)((number_of_states - 1.0) + 1.0) - 1) {
    d0 = rt_powd_snf(2.0, k) - 1.0;
    l = 0;
    for (b_l = 0; b_l < (int)(d0 + 1.0); b_l++) {
      l = b_l;
      nxt_stat(n, b_l, L, k, &next_state, dec_output_bin);
      input->data[n + input->size[0] * ((int)(next_state + 1.0) - 1)] = b_l;

      /* 输入矩阵，由当前状态和下一状态确定输入的十进制表示 */
      for (i0 = 0; i0 < 2; i0++) {
        for (b_k = 0; b_k < 3; b_k++) {
          b[b_k + 3 * i0] = G[i0 + (b_k << 1)];
        }
      }

      if (dec_output_bin->size[1] == 1) {
        for (i0 = 0; i0 < 2; i0++) {
          y[i0] = 0.0;
          for (b_k = 0; b_k < 3; b_k++) {
            mtmp = y[i0] + dec_output_bin->data[b_k] * b[b_k + 3 * i0];
            y[i0] = mtmp;
          }
        }
      } else {
        b_k = dec_output_bin->size[1];
        for (i0 = 0; i0 < 2; i0++) {
          y[i0] = 0.0;
        }

        for (nx = 0; nx < 2; nx++) {
          for (ic = nx; ic + 1 <= nx + 1; ic++) {
            y[ic] = 0.0;
          }
        }

        br = 0;
        for (nx = 0; nx < 2; nx++) {
          ar = -1;
          i0 = br + b_k;
          for (ib = br; ib + 1 <= i0; ib++) {
            if (b[ib] != 0.0) {
              idx = ar;
              for (ic = nx; ic + 1 <= nx + 1; ic++) {
                idx++;
                y[ic] += b[ib] * dec_output_bin->data[idx];
              }
            }

            ar++;
          }

          br += b_k;
        }
      }

      nextstate->data[n + nextstate->size[0] * b_l] = next_state;

      /* 转移矩阵，下一状态的十进制表示 */
      /*  将二进制数转化为十进制数 */
      b_rem(y, binary_output);
      mtmp = 0.0;
      for (i0 = 0; i0 < 2; i0++) {
        mtmp += binary_output[i0] * (2.0 - (double)i0);
      }

      output->data[n + output->size[0] * b_l] = mtmp;

      /* 输出矩阵，输出的十进制表示 */
    }

    n++;
  }

  emxInit_real_T(&state_metric, 2);

  /* state_metric状态距离矩阵（状态数x2）,第二列表明到达下一状态的最小距离，第一列表明当前状态的最小距离. */
  i0 = state_metric->size[0] * state_metric->size[1];
  state_metric->size[0] = (int)number_of_states;
  state_metric->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)state_metric, i0, (int)sizeof(double));
  ic = (int)number_of_states << 1;
  for (i0 = 0; i0 < ic; i0++) {
    state_metric->data[i0] = 0.0;
  }

  emxInit_real_T(&survivor_state, 2);

  /* 记录当前级各节点的最小欧式距离和下一级各节点的欧式距离 */
  /* 级数 */
  /* 矩阵的每一行代表编码器的一个输出 */
  /* survivor_state(状态数x结点深度)幸存状态矩阵：（i,j）表示第j个时刻时第i个状态是由谁过来的。 */
  i0 = survivor_state->size[0] * survivor_state->size[1];
  survivor_state->size[0] = (int)number_of_states;
  survivor_state->size[1] = 1003;
  emxEnsureCapacity((emxArray__common *)survivor_state, i0, (int)sizeof(double));
  ic = (int)number_of_states * 1003;
  for (i0 = 0; i0 < ic; i0++) {
    survivor_state->data[i0] = 0.0;
  }

  emxInit_real_T(&survivor_first_bit, 2);

  /* 通过网格的最佳路径的矩阵 */
  i0 = survivor_first_bit->size[0] * survivor_first_bit->size[1];
  survivor_first_bit->size[0] = (int)number_of_states;
  survivor_first_bit->size[1] = 1003;
  emxEnsureCapacity((emxArray__common *)survivor_first_bit, i0, (int)sizeof
                    (double));
  ic = (int)number_of_states * 1003;
  for (i0 = 0; i0 < ic; i0++) {
    survivor_first_bit->data[i0] = 0.0;
  }

  emxInit_real_T(&survivor_second_bit, 2);
  i0 = survivor_second_bit->size[0] * survivor_second_bit->size[1];
  survivor_second_bit->size[0] = (int)number_of_states;
  survivor_second_bit->size[1] = 1003;
  emxEnsureCapacity((emxArray__common *)survivor_second_bit, i0, (int)sizeof
                    (double));
  ic = (int)number_of_states * 1003;
  for (i0 = 0; i0 < ic; i0++) {
    survivor_second_bit->data[i0] = 0.0;
  }

  /* 开始无尾信道输出的解码（即去掉结尾(L-1)个零） */
  /* body encoder */
  m = 0;
  emxInit_real_T(&b_state_metric, 2);
  while (m <= (int)((1002.0 - L) + 1.0) - 1) {
    i0 = dec_output_bin->size[0] * dec_output_bin->size[1];
    dec_output_bin->size[0] = 1;
    dec_output_bin->size[1] = (int)number_of_states;
    emxEnsureCapacity((emxArray__common *)dec_output_bin, i0, (int)sizeof(double));
    ic = (int)number_of_states;
    for (i0 = 0; i0 < ic; i0++) {
      dec_output_bin->data[i0] = 0.0;
    }

    /* 标志各状态有否被访问过 */
    /* flag(1x状态数)：第j位表示在一个时刻中第j个状态是否到达过，1表示到达过。 */
    /* 确定状态转移的步长，起始阶段的步长与中间态的步长不同 */
    if (1.0 + (double)m <= L) {
      step = rt_powd_snf(2.0, (L - (1.0 + (double)m)) * k);

      /* 表示第i级各状态的最小差值 */
    } else {
      step = 1.0;
    }

    /* 进行加－比－选 */
    i0 = (int)(((number_of_states - 1.0) + step) / step);
    for (n = 0; n < i0; n++) {
      b_n = (double)n * step;
      d0 = rt_powd_snf(2.0, k) - 1.0;
      l = 0;
      for (b_l = 0; b_l < (int)(d0 + 1.0); b_l++) {
        l = b_l;
        b_deci2bin(output->data[((int)(b_n + 1.0) + output->size[0] * b_l) - 1],
                   binary_output);

        /* 将理想输出转化为二进制 */
        parallel_binary_output[0] = 0.0;
        parallel_binary_output[1] = 0.0;

        /* 并行转移路径 */
        parallel_binary_output_two[0] = 0.0;
        parallel_binary_output_two[1] = 1.0;
        parallel_binary_output_three[0] = 1.0;
        parallel_binary_output_three[1] = 0.0;

        /* 并行转移路径 */
        parallel_binary_output_four[0] = 1.0;
        parallel_binary_output_four[1] = 1.0;
        for (b_k = 0; b_k < 2; b_k++) {
          parallel_binary_output[b_k + 2] = binary_output[b_k];
          parallel_binary_output_two[b_k + 2] = binary_output[b_k];
          parallel_binary_output_three[b_k + 2] = binary_output[b_k];
          parallel_binary_output_four[b_k + 2] = binary_output[b_k];
        }

        /*  将二进制数转化为十进制数 */
        mtmp = 0.0;
        for (b_k = 0; b_k < 4; b_k++) {
          mtmp += parallel_binary_output[b_k] * (double)b_b[b_k];
        }

        /*  将二进制数转化为十进制数 */
        c_a = 0.0;
        for (b_k = 0; b_k < 4; b_k++) {
          c_a += parallel_binary_output_two[b_k] * (double)b_b[b_k];
        }

        /*  将二进制数转化为十进制数 */
        d_a = 0.0;
        for (b_k = 0; b_k < 4; b_k++) {
          d_a += parallel_binary_output_three[b_k] * (double)b_b[b_k];
        }

        /*  将二进制数转化为十进制数 */
        d1 = 0.0;
        for (b_k = 0; b_k < 4; b_k++) {
          d1 += parallel_binary_output_four[b_k] * (double)b_b[b_k];
        }

        /*  if x==y */
        /*      distance=0; */
        /*  else  */
        /*      distance=1; */
        /*  end */
        /*  欧式距离 */
        a = channel_output[m].re - s16[(int)(mtmp + 1.0) - 1].re;
        b_a = channel_output[m].im - s16[(int)(mtmp + 1.0) - 1].im;

        /*  switch (y) */
        /*  case 0 */
        /*      if x==0 */
        /*          distance=0.0458; */
        /*      end */
        /*      if x==1 */
        /*          distance=2; */
        /*      end */
        /*      if x==2; */
        /*          distance==1.0458; */
        /*      end */
        /*  case 1 */
        /*      if x==0 */
        /*          distance==2; */
        /*      end */
        /*      if x==1 */
        /*          distance=0.0458; */
        /*      end */
        /*      if x==2 */
        /*          distance=1.0458; */
        /*      end */
        /*  otherwise, */
        /*      break; */
        /*  end */
        /*  if x==y */
        /*      distance=0; */
        /*  else  */
        /*      distance=1; */
        /*  end */
        /*  欧式距离 */
        e_a = channel_output[m].re - s16[(int)(c_a + 1.0) - 1].re;
        f_a = channel_output[m].im - s16[(int)(c_a + 1.0) - 1].im;

        /*  switch (y) */
        /*  case 0 */
        /*      if x==0 */
        /*          distance=0.0458; */
        /*      end */
        /*      if x==1 */
        /*          distance=2; */
        /*      end */
        /*      if x==2; */
        /*          distance==1.0458; */
        /*      end */
        /*  case 1 */
        /*      if x==0 */
        /*          distance==2; */
        /*      end */
        /*      if x==1 */
        /*          distance=0.0458; */
        /*      end */
        /*      if x==2 */
        /*          distance=1.0458; */
        /*      end */
        /*  otherwise, */
        /*      break; */
        /*  end */
        /*  if x==y */
        /*      distance=0; */
        /*  else  */
        /*      distance=1; */
        /*  end */
        /*  欧式距离 */
        g_a = channel_output[m].re - s16[(int)(d_a + 1.0) - 1].re;
        c_a = channel_output[m].im - s16[(int)(d_a + 1.0) - 1].im;

        /*  switch (y) */
        /*  case 0 */
        /*      if x==0 */
        /*          distance=0.0458; */
        /*      end */
        /*      if x==1 */
        /*          distance=2; */
        /*      end */
        /*      if x==2; */
        /*          distance==1.0458; */
        /*      end */
        /*  case 1 */
        /*      if x==0 */
        /*          distance==2; */
        /*      end */
        /*      if x==1 */
        /*          distance=0.0458; */
        /*      end */
        /*      if x==2 */
        /*          distance=1.0458; */
        /*      end */
        /*  otherwise, */
        /*      break; */
        /*  end */
        /*  if x==y */
        /*      distance=0; */
        /*  else  */
        /*      distance=1; */
        /*  end */
        /*  欧式距离 */
        d_a = channel_output[m].re - s16[(int)(d1 + 1.0) - 1].re;
        mtmp = channel_output[m].im - s16[(int)(d1 + 1.0) - 1].im;

        /*  switch (y) */
        /*  case 0 */
        /*      if x==0 */
        /*          distance=0.0458; */
        /*      end */
        /*      if x==1 */
        /*          distance=2; */
        /*      end */
        /*      if x==2; */
        /*          distance==1.0458; */
        /*      end */
        /*  case 1 */
        /*      if x==0 */
        /*          distance==2; */
        /*      end */
        /*      if x==1 */
        /*          distance=0.0458; */
        /*      end */
        /*      if x==2 */
        /*          distance=1.0458; */
        /*      end */
        /*  otherwise, */
        /*      break; */
        /*  end */
        parallel_diatance_output[0] = a * a + b_a * b_a;
        parallel_diatance_output[1] = e_a * e_a + f_a * f_a;
        parallel_diatance_output[2] = g_a * g_a + c_a * c_a;
        parallel_diatance_output[3] = d_a * d_a + mtmp * mtmp;
        nx = 1;
        mtmp = parallel_diatance_output[0];
        if (rtIsNaN(parallel_diatance_output[0])) {
          idx = 2;
          exitg2 = false;
          while ((!exitg2) && (idx < 5)) {
            nx = idx;
            if (!rtIsNaN(parallel_diatance_output[idx - 1])) {
              mtmp = parallel_diatance_output[idx - 1];
              exitg2 = true;
            } else {
              idx++;
            }
          }
        }

        if (nx < 4) {
          while (nx + 1 < 5) {
            if (parallel_diatance_output[nx] < mtmp) {
              mtmp = parallel_diatance_output[nx];
            }

            nx++;
          }
        }

        /* 计算欧式距离 */
        idx = 0;
        for (b_k = 0; b_k < 2; b_k++) {
          ii_size[b_k] = iv0[b_k];
        }

        nx = 1;
        exitg1 = false;
        while ((!exitg1) && (nx < 5)) {
          guard1 = false;
          if (parallel_diatance_output[nx - 1] == mtmp) {
            idx++;
            ii_data[idx - 1] = (signed char)nx;
            if (idx >= 4) {
              exitg1 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            nx++;
          }
        }

        if (1 > idx) {
          ic = 0;
        } else {
          ic = idx;
        }

        if (1 > idx) {
          ii_size[1] = 0;
        } else {
          ii_size[1] = idx;
        }

        for (b_k = 0; b_k < 4; b_k++) {
          b_parallel_binary_output[b_k << 2] = parallel_binary_output[b_k];
          b_parallel_binary_output[1 + (b_k << 2)] =
            parallel_binary_output_two[b_k];
          b_parallel_binary_output[2 + (b_k << 2)] =
            parallel_binary_output_three[b_k];
          b_parallel_binary_output[3 + (b_k << 2)] =
            parallel_binary_output_four[b_k];
        }

        for (b_k = 0; b_k < 4; b_k++) {
          for (nx = 0; nx < ic; nx++) {
            parallel_double_bit_output_data[nx + ii_size[1] * b_k] =
              b_parallel_binary_output[(ii_data[ii_size[0] * nx] + (b_k << 2)) -
              1];
          }
        }

        /* 在AWGN信道下，最大似然估计转化为求最小欧式距离 */
        /* 如果下一状态度量距离大于当前距离加欧式距离，或是下一状态未被遍历过则设为当前状态下一状态的幸存状态，当前距离加欧式距离设为下一状态的距离   */
        if ((state_metric->data[((int)(nextstate->data[((int)(b_n + 1.0) +
                nextstate->size[0] * b_l) - 1] + 1.0) + state_metric->size[0]) -
             1] > state_metric->data[(int)(b_n + 1.0) - 1] + mtmp) || ((int)
             dec_output_bin->data[(int)(nextstate->data[((int)(b_n + 1.0) +
               nextstate->size[0] * b_l) - 1] + 1.0) - 1] == 0)) {
          state_metric->data[((int)(nextstate->data[((int)(b_n + 1.0) +
            nextstate->size[0] * b_l) - 1] + 1.0) + state_metric->size[0]) - 1] =
            state_metric->data[(int)(b_n + 1.0) - 1] + mtmp;

          /* 更改欧式距离 */
          survivor_state->data[((int)(nextstate->data[((int)(b_n + 1.0) +
            nextstate->size[0] * b_l) - 1] + 1.0) + survivor_state->size[0] *
                                ((int)((1.0 + (double)m) + 1.0) - 1)) - 1] = b_n;

          /* 更改幸存路径 */
          survivor_first_bit->data[((int)(nextstate->data[((int)(b_n + 1.0) +
            nextstate->size[0] * b_l) - 1] + 1.0) + survivor_first_bit->size[0] *
            ((int)((1.0 + (double)m) + 1.0) - 1)) - 1] =
            parallel_double_bit_output_data[0];
          survivor_second_bit->data[((int)(nextstate->data[((int)(b_n + 1.0) +
            nextstate->size[0] * b_l) - 1] + 1.0) + survivor_second_bit->size[0]
            * ((int)((1.0 + (double)m) + 1.0) - 1)) - 1] =
            parallel_double_bit_output_data[1];
          dec_output_bin->data[(int)(nextstate->data[((int)(b_n + 1.0) +
            nextstate->size[0] * b_l) - 1] + 1.0) - 1] = 1.0;
        }
      }
    }

    nx = state_metric->size[0];
    i0 = b_state_metric->size[0] * b_state_metric->size[1];
    b_state_metric->size[0] = nx;
    b_state_metric->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)b_state_metric, i0, (int)sizeof(double));
    for (i0 = 0; i0 < 2; i0++) {
      for (b_k = 0; b_k < nx; b_k++) {
        b_state_metric->data[b_k + b_state_metric->size[0] * i0] =
          state_metric->data[b_k + state_metric->size[0] * (1 - i0)];
      }
    }

    i0 = state_metric->size[0] * state_metric->size[1];
    state_metric->size[0] = b_state_metric->size[0];
    state_metric->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)state_metric, i0, (int)sizeof(double));
    for (i0 = 0; i0 < 2; i0++) {
      ic = b_state_metric->size[0];
      for (b_k = 0; b_k < ic; b_k++) {
        state_metric->data[b_k + state_metric->size[0] * i0] =
          b_state_metric->data[b_k + b_state_metric->size[0] * i0];
      }
    }

    /* 将两列互换，目的是上一次计算的下一状态距离变为下一次计算的上一状态距离 */
    m++;
  }

  emxFree_real_T(&b_state_metric);

  /* 尾部信道输出的解码，最后的那些级没有包含所有的状态 */
  /* 开始尾部信道输出解码 */
  /* 归零计算 */
  /* tailer encoder */
  m = 0;
  emxInit_real_T(&c_state_metric, 2);
  while (m <= (int)(1002.0 + (1.0 - ((1002.0 - L) + 2.0))) - 1) {
    b_m = ((1002.0 - L) + 2.0) + (double)m;
    i0 = dec_output_bin->size[0] * dec_output_bin->size[1];
    dec_output_bin->size[0] = 1;
    dec_output_bin->size[1] = (int)number_of_states;
    emxEnsureCapacity((emxArray__common *)dec_output_bin, i0, (int)sizeof(double));
    ic = (int)number_of_states;
    for (i0 = 0; i0 < ic; i0++) {
      dec_output_bin->data[i0] = 0.0;
    }

    /* flag(1x状态数)：第j位表示在一个时刻中第j个状态是否到达过，1表示到达过。 */
    /* 所剩余的状态数 */
    /* 认为尾部的数据为零，每次状态转移，只有输入为零一条路径     */
    d0 = number_of_states / rt_powd_snf(2.0, (((b_m - 1002.0) + L) - 2.0) * k) -
      1.0;
    for (n = 0; n < (int)(d0 + 1.0); n++) {
      b_deci2bin(output->data[n], binary_output);
      parallel_binary_output[0] = 0.0;
      parallel_binary_output[1] = 0.0;
      for (i0 = 0; i0 < 2; i0++) {
        parallel_binary_output[i0 + 2] = binary_output[i0];
      }

      /* 并行转移路径 */
      /*  将二进制数转化为十进制数 */
      mtmp = 0.0;
      for (i0 = 0; i0 < 4; i0++) {
        mtmp += parallel_binary_output[i0] * (double)b_b[i0];
      }

     
      a = channel_output[(int)b_m - 1].re - s16[(int)(mtmp + 1.0) - 1].re;
      b_a = channel_output[(int)b_m - 1].im - s16[(int)(mtmp + 1.0) - 1].im;
      branch_metric = a * a + b_a * b_a;

     
      if ((state_metric->data[((int)(nextstate->data[n] + 1.0) +
            state_metric->size[0]) - 1] > state_metric->data[n] + branch_metric)
          || ((int)dec_output_bin->data[(int)(nextstate->data[n] + 1.0) - 1] ==
              0)) {
        state_metric->data[((int)(nextstate->data[n] + 1.0) + state_metric->
                            size[0]) - 1] = state_metric->data[n] +
          branch_metric;

        /* 在第2列保存下一状态的距离 */
        survivor_state->data[((int)(nextstate->data[n] + 1.0) +
                              survivor_state->size[0] * (int)b_m) - 1] = n;

        /* 更改幸存路径 */
        survivor_first_bit->data[((int)(nextstate->data[n + nextstate->size[0] *
          l] + 1.0) + survivor_first_bit->size[0] * (int)b_m) - 1] = 0.0;
        survivor_second_bit->data[((int)(nextstate->data[n + nextstate->size[0] *
          l] + 1.0) + survivor_second_bit->size[0] * (int)b_m) - 1] = 0.0;
        dec_output_bin->data[(int)(nextstate->data[n] + 1.0) - 1] = 1.0;
      }
    }

    nx = state_metric->size[0];
    i0 = c_state_metric->size[0] * c_state_metric->size[1];
    c_state_metric->size[0] = nx;
    c_state_metric->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)c_state_metric, i0, (int)sizeof(double));
    for (i0 = 0; i0 < 2; i0++) {
      for (b_k = 0; b_k < nx; b_k++) {
        c_state_metric->data[b_k + c_state_metric->size[0] * i0] =
          state_metric->data[b_k + state_metric->size[0] * (1 - i0)];
      }
    }

    i0 = state_metric->size[0] * state_metric->size[1];
    state_metric->size[0] = c_state_metric->size[0];
    state_metric->size[1] = 2;
    emxEnsureCapacity((emxArray__common *)state_metric, i0, (int)sizeof(double));
    for (i0 = 0; i0 < 2; i0++) {
      ic = c_state_metric->size[0];
      for (b_k = 0; b_k < ic; b_k++) {
        state_metric->data[b_k + state_metric->size[0] * i0] =
          c_state_metric->data[b_k + c_state_metric->size[0] * i0];
      }
    }

    /* 将两列互换，目的是上一次计算的下一状态距离变为下一次计算的上一状态距离 */
    m++;
  }

  emxFree_real_T(&c_state_metric);
  emxFree_real_T(&state_metric);
  emxFree_real_T(&output);
  emxFree_real_T(&nextstate);

  /* 开始回溯最佳路径，从最佳路径中找出解码 */
  /* state_sequence(1x结点深度)矩阵，1～dep分别记载各个阶段的路径（即前一状态数）。1～dep分别记载各个阶段的路径（即前一状态数）。 */
  /* 从最佳路径中产生解码 */
  /* 由后到前得到各级的状态 */
  memset(&state_sequence[0], 0, 1003U * sizeof(double));
  state_sequence[1001] = survivor_state->data[survivor_state->size[0] * 1002];

  /* 开始回溯最佳路径 */
  for (m = 0; m < 1002; m++) {
    state_sequence[1001 - m] = survivor_state->data[((int)(state_sequence[1002 -
      m] + 1.0) + survivor_state->size[0] * (1002 - m)) - 1];
    first_bit_sequence[1001 - m] = survivor_first_bit->data[((int)
      (state_sequence[1002 - m] + 1.0) + survivor_first_bit->size[0] * (1002 - m))
      - 1];
    second_bit_sequence[1001 - m] = survivor_second_bit->data[((int)
      (state_sequence[1002 - m] + 1.0) + survivor_second_bit->size[0] * (1002 -
      m)) - 1];
  }

  emxFree_real_T(&survivor_second_bit);
  emxFree_real_T(&survivor_first_bit);
  emxFree_real_T(&survivor_state);
  emxInit_real_T(&decoder_output_matrix, 2);

  /* 解码输出矩阵，输出码元个数为结点深度减L */
  i0 = decoder_output_matrix->size[0] * decoder_output_matrix->size[1];
  decoder_output_matrix->size[0] = (int)k;
  decoder_output_matrix->size[1] = (int)((1002.0 - L) + 1.0);
  emxEnsureCapacity((emxArray__common *)decoder_output_matrix, i0, (int)sizeof
                    (double));
  ic = (int)k * (int)((1002.0 - L) + 1.0);
  for (i0 = 0; i0 < ic; i0++) {
    decoder_output_matrix->data[i0] = 0.0;
  }

  /* 由各级的状态和输入矩阵得到各级的输入 */
  for (m = 0; m < (int)((1002.0 - L) + 1.0); m++) {
    /* 解出码元 */
    deci2bin(input->data[((int)(state_sequence[m] + 1.0) + input->size[0] *
                          ((int)(state_sequence[(int)((1.0 + (double)m) + 1.0) -
                1] + 1.0) - 1)) - 1], k, dec_output_bin);

    /* 按k将输出码元转换 */
    if (1.0 > k) {
      i0 = 1;
      b_k = 1;
      nx = 0;
    } else {
      i0 = (int)k;
      b_k = -1;
      nx = 1;
    }

    ic = div_s32_floor(nx - i0, b_k);
    for (nx = 0; nx <= ic; nx++) {
      decoder_output_matrix->data[nx + decoder_output_matrix->size[0] * m] =
        dec_output_bin->data[(i0 + b_k * nx) - 1];
    }
  }

  emxFree_real_T(&dec_output_bin);
  if (1.0 > (1002.0 - L) + 1.0) {
    ic = 0;
  } else {
    ic = (int)((1002.0 - L) + 1.0);
  }

  if (1.0 > (1002.0 - L) + 1.0) {
    nx = 0;
  } else {
    nx = (int)((1002.0 - L) + 1.0);
  }

  /* 由各级的输入的得到输入序列 */
  if (!(ic == 0)) {
    idx = ic;
  } else if (!(nx == 0)) {
    idx = nx;
  } else if (!((decoder_output_matrix->size[0] == 0) ||
               (decoder_output_matrix->size[1] == 0))) {
    idx = decoder_output_matrix->size[1];
  } else {
    idx = 0;
    if (decoder_output_matrix->size[1] > 0) {
      idx = decoder_output_matrix->size[1];
    }
  }

  empty_non_axis_sizes = (idx == 0);
  if (empty_non_axis_sizes || (!(ic == 0))) {
    br = 1;
  } else {
    br = 0;
  }

  if (empty_non_axis_sizes || (!(nx == 0))) {
    ar = 1;
  } else {
    ar = 0;
  }

  if (empty_non_axis_sizes || (!((decoder_output_matrix->size[0] == 0) ||
        (decoder_output_matrix->size[1] == 0)))) {
    ib = decoder_output_matrix->size[0];
  } else {
    ib = 0;
  }

  for (i0 = 0; i0 < ic; i0++) {
    first_bit_sequence_data[i0] = first_bit_sequence[i0];
  }

  for (i0 = 0; i0 < nx; i0++) {
    second_bit_sequence_data[i0] = second_bit_sequence[i0];
  }

  i0 = input->size[0] * input->size[1];
  input->size[0] = (br + ar) + ib;
  input->size[1] = idx;
  emxEnsureCapacity((emxArray__common *)input, i0, (int)sizeof(double));
  for (i0 = 0; i0 < idx; i0++) {
    for (b_k = 0; b_k < br; b_k++) {
      input->data[b_k + input->size[0] * i0] = first_bit_sequence_data[b_k + br *
        i0];
    }
  }

  for (i0 = 0; i0 < idx; i0++) {
    for (b_k = 0; b_k < ar; b_k++) {
      input->data[(b_k + br) + input->size[0] * i0] =
        second_bit_sequence_data[b_k + ar * i0];
    }
  }

  for (i0 = 0; i0 < idx; i0++) {
    for (b_k = 0; b_k < ib; b_k++) {
      input->data[((b_k + br) + ar) + input->size[0] * i0] =
        decoder_output_matrix->data[b_k + ib * i0];
    }
  }

  emxFree_real_T(&decoder_output_matrix);
  nx = input->size[0] * input->size[1];
  i0 = decoder_output->size[0] * decoder_output->size[1];
  decoder_output->size[0] = 1;
  decoder_output->size[1] = (int)(3.0 * k * ((1002.0 - L) + 1.0));
  emxEnsureCapacity((emxArray__common *)decoder_output, i0, (int)sizeof(double));
  for (b_k = 0; b_k + 1 <= nx; b_k++) {
    decoder_output->data[b_k] = input->data[b_k];
  }

  emxFree_real_T(&input);
}

/*
 * File trailer for vi_de.c
 *
 * [EOF]
 */
