/*
 * File: main.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 15-Jul-2016 10:32:49
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "vi_de.h"
#include "main.h"
#include "vi_de_terminate.h"
#include "vi_de_emxAPI.h"
#include "vi_de_initialize.h"
#include <stdio.h>

/* Function Declarations */
static void argInit_1x1002_creal_T(creal_T result[1002]);
static void argInit_1x16_creal_T(creal_T result[16]);
static void argInit_2x3_real_T(double result[6]);
static creal_T argInit_creal_T(void);
static double argInit_real_T(void);
static void main_vi_de(void);

/* Function Definitions */

/*
 * Arguments    : creal_T result[1002]
 * Return Type  : void
 */
static void argInit_1x1002_creal_T(creal_T result[1002])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  //for (idx1 = 0; idx1 < 1002; idx1++) {
  //  /* Set the value of the array element.
  //     Change this value to the value that the application requires. */
  //  result[idx1] = argInit_creal_T();
  //}
  FILE *p;
  errno_t  fp = fopen_s(&p, "result.txt", "r");
  if (!p)
  {
	  printf("can't open file\n");
	  return ;
  }
  for (idx1 = 0; idx1 < 1002; idx1++)
  {
	  int d = 3;
	  fscanf_s(p, "%d", &d);
	  result[idx1].re = d;
	  fscanf_s(p, "%d", &d);
	  result[idx1].im = d;
  }
  fclose(p);
}

/*
 * Arguments    : creal_T result[16]
 * Return Type  : void
 */
static void argInit_1x16_creal_T(creal_T s16[16])
{
  //int idx1;

  ///* Loop over the array to initialize each element. */
  //for (idx1 = 0; idx1 < 16; idx1++) {
  //  /* Set the value of the array element.
  //     Change this value to the value that the application requires. */
  //  result[idx1] = argInit_creal_T();
  //}
  s16[0].re = -3;
  s16[0].im = 3;
  s16[4].re = 1;
  s16[4].im = -1;
  s16[8].re = -3;
  s16[8].im = -1;
  s16[12].re = 1;
  s16[12].im = 3;
  s16[1].re = -1;
  s16[1].im = 1;
  s16[5].re = 3;
  s16[5].im = -3;
  s16[9].re = -1;
  s16[9].im = -3;
  s16[13].re = 3;
  s16[13].im = 1;
  s16[2].re = -1;
  s16[2].im = -1;
  s16[6].re = 3;
  s16[6].im = 3;
  s16[10].re = -1;
  s16[10].im = 3;
  s16[14].re = 3;
  s16[14].im = -1;
  s16[3].re = -3;
  s16[3].im = -3;
  s16[7].re = 1;
  s16[7].im = 1;
  s16[11].re = -3;
  s16[11].im = 1;
  s16[15].re = 1;
  s16[15].im = -3;

}

/*
 * Arguments    : double result[6]
 * Return Type  : void
 */
static void argInit_2x3_real_T(double result[6])
{
  int idx0;
  int idx1;
  int aa[6] = { 1, 0, 0, 1, 1, 0 };
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    for (idx1 = 0; idx1 < 3; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
		result[idx0 + (idx1 << 1)] = aa[idx0 + (idx1 << 1)];
		  //argInit_real_T();
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : creal_T
 */
static creal_T argInit_creal_T(void)
{
  creal_T result;

  /* Set the value of the complex variable.
     Change this value to the value that the application requires. */
  result.re = argInit_real_T();
  result.im = argInit_real_T();
  return result;
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_vi_de(void)
{
  emxArray_real_T *decoder_output;
  double dv0[6];
  creal_T dcv0[1002];
  creal_T dcv1[16];
  emxInitArray_real_T(&decoder_output, 2);

  /* Initialize function 'vi_de' input arguments. */
  /* Initialize function input argument 'G'. */
  /* Initialize function input argument 'channel_output'. */
  /* Initialize function input argument 's16'. */
  /* Call the entry-point 'vi_de'. */
  argInit_2x3_real_T(dv0);
  argInit_1x1002_creal_T(dcv0);
  argInit_1x16_creal_T(dcv1);
  double k_code = 1.0;
  vi_de(dv0, k_code, dcv0, dcv1, decoder_output);
  FILE * p;
  errno_t  fp = fopen_s(&p, "decode.txt", "w");
  if (!p)
  {
	  printf("can't open file\n");
	  return;
  }
  for (int idx = 0; idx < 3000; idx++)
  {
	  fprintf(p, "%d", (int)decoder_output->data[idx]);
	  putc('\n', p);
  }
  emxDestroyArray_real_T(decoder_output);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  vi_de_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_vi_de();

  /* Terminate the application.
     You do not need to do this more than one time. */
  vi_de_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
