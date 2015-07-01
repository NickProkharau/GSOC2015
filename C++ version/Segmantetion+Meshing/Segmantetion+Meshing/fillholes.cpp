
// Include files
#include "stdafx.h"
#include "rt_nonfinite.h"
#include "cgalv2m.h"
#include "corner_points.h"
#include "fillholes.h"
#include "fillimg.h"
#include "meshing.h"
#include "read_mri_simple.h"
#include "readmedit.h"
#include "rot_matrix.h"
#include "saveinr.h"
#include "segmentation.h"
#include "sphere.h"
#include "volumesegment.h"
#include "warp_apply.h"
#include "Segmentation1_emxutil.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : const emxArray_real_T *input
//                emxArray_boolean_T *output
// Return Type  : void
//
void b_fillholes(const emxArray_real_T *input, emxArray_boolean_T *output)
{
  emxArray_real_T *output1;
  int i31;
  int loop_ub;
  unsigned int dim[7];
  int i;
  emxArray_real_T *slice1;
  emxArray_real_T *r12;
  int k;
  boolean_T exitg2;
  int i32;
  int slice1_idx_0;
  int input_idx_1;
  int sqsz[3];
  int j;
  emxArray_real_T *output2;
  boolean_T exitg1;
  emxArray_real_T *output3;
  emxArray_real_T *b_slice1;
  int b_loop_ub;
  d_emxInit_real_T(&output1, 7);
  i31 = output1->size[0] * output1->size[1] * output1->size[2] * output1->size[3]
    * output1->size[4] * output1->size[5] * output1->size[6];
  output1->size[0] = input->size[0];
  output1->size[1] = input->size[1];
  output1->size[2] = input->size[2];
  output1->size[3] = input->size[3];
  output1->size[4] = input->size[4];
  output1->size[5] = input->size[5];
  output1->size[6] = input->size[6];
  emxEnsureCapacity((emxArray__common *)output1, i31, (int)sizeof(double));
  loop_ub = input->size[0] * input->size[1] * input->size[2] * input->size[3] *
    input->size[4] * input->size[5] * input->size[6];
  for (i31 = 0; i31 < loop_ub; i31++) {
    output1->data[i31] = input->data[i31];
  }

  for (i31 = 0; i31 < 7; i31++) {
    dim[i31] = (unsigned int)input->size[i31];
  }

  i = 0;
  b_emxInit_real_T(&slice1, 2);
  b_emxInit_real_T(&r12, 2);
  while (i <= (int)dim[0] - 1) {
    k = 3;
    exitg2 = false;
    while ((!exitg2) && (k > 2)) {
      i31 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
        input->size[6];
      if (i31 == 1) {
        k = 2;
      } else {
        exitg2 = true;
      }
    }

    if (k <= 2) {
      i31 = input->size[1];
      i32 = slice1->size[0] * slice1->size[1];
      slice1->size[0] = 1;
      slice1->size[1] = i31;
      emxEnsureCapacity((emxArray__common *)slice1, i32, (int)sizeof(double));
      i31 = input->size[1];
      i32 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
        input->size[6];
      i31 *= i32;
      for (k = 0; k + 1 <= i31; k++) {
        i32 = input->size[1];
        slice1_idx_0 = input->size[0];
        input_idx_1 = input->size[1];
        slice1->data[k] = input->data[(i + slice1_idx_0 * (k % i32)) +
          slice1_idx_0 * input_idx_1 * (k / i32)];
      }
    } else {
      for (i31 = 0; i31 < 3; i31++) {
        sqsz[i31] = 1;
      }

      j = 0;
      i31 = input->size[1];
      if (i31 != 1) {
        i31 = input->size[1];
        sqsz[0] = i31;
        j = 1;
      }

      i31 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
        input->size[6];
      if (i31 != 1) {
        i31 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
          input->size[6];
        sqsz[j] = i31;
      }

      i31 = slice1->size[0] * slice1->size[1];
      slice1->size[0] = sqsz[0];
      slice1->size[1] = sqsz[1];
      emxEnsureCapacity((emxArray__common *)slice1, i31, (int)sizeof(double));
      i31 = input->size[1];
      i32 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
        input->size[6];
      i31 *= i32;
      for (k = 0; k + 1 <= i31; k++) {
        i32 = input->size[1];
        slice1_idx_0 = input->size[0];
        input_idx_1 = input->size[1];
        slice1->data[k] = input->data[(i + slice1_idx_0 * (k % i32)) +
          slice1_idx_0 * input_idx_1 * (k / i32)];
      }
    }

    e_fillimg(slice1, r12);
    k = output1->size[0];
    slice1_idx_0 = r12->size[0];
    j = r12->size[1];
    loop_ub = j - 1;
    for (i31 = 0; i31 <= loop_ub; i31++) {
      output1->data[i + k * i31] = r12->data[i + slice1_idx_0 * i31];
    }

    i++;
  }

  d_emxInit_real_T(&output2, 7);
  i31 = output2->size[0] * output2->size[1] * output2->size[2] * output2->size[3]
    * output2->size[4] * output2->size[5] * output2->size[6];
  output2->size[0] = input->size[0];
  output2->size[1] = input->size[1];
  output2->size[2] = input->size[2];
  output2->size[3] = input->size[3];
  output2->size[4] = input->size[4];
  output2->size[5] = input->size[5];
  output2->size[6] = input->size[6];
  emxEnsureCapacity((emxArray__common *)output2, i31, (int)sizeof(double));
  loop_ub = input->size[0] * input->size[1] * input->size[2] * input->size[3] *
    input->size[4] * input->size[5] * input->size[6];
  for (i31 = 0; i31 < loop_ub; i31++) {
    output2->data[i31] = input->data[i31];
  }

  for (i = 0; i < (int)dim[1]; i++) {
    k = 3;
    exitg1 = false;
    while ((!exitg1) && (k > 2)) {
      i31 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
        input->size[6];
      if (i31 == 1) {
        k = 2;
      } else {
        exitg1 = true;
      }
    }

    if (k <= 2) {
      i31 = input->size[0];
      i32 = slice1->size[0] * slice1->size[1];
      slice1->size[0] = i31;
      slice1->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)slice1, i32, (int)sizeof(double));
      i31 = input->size[0];
      i32 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
        input->size[6];
      i31 *= i32;
      for (k = 0; k + 1 <= i31; k++) {
        i32 = input->size[0];
        slice1_idx_0 = input->size[0];
        input_idx_1 = input->size[1];
        slice1->data[k] = input->data[(k % i32 + slice1_idx_0 * i) +
          slice1_idx_0 * input_idx_1 * (k / i32)];
      }
    } else {
      for (i31 = 0; i31 < 3; i31++) {
        sqsz[i31] = 1;
      }

      j = 0;
      i31 = input->size[0];
      if (i31 != 1) {
        i31 = input->size[0];
        sqsz[0] = i31;
        j = 1;
      }

      i31 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
        input->size[6];
      if (i31 != 1) {
        i31 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
          input->size[6];
        sqsz[j] = i31;
      }

      i31 = slice1->size[0] * slice1->size[1];
      slice1->size[0] = sqsz[0];
      slice1->size[1] = sqsz[1];
      emxEnsureCapacity((emxArray__common *)slice1, i31, (int)sizeof(double));
      i31 = input->size[0];
      i32 = input->size[2] * input->size[3] * input->size[4] * input->size[5] *
        input->size[6];
      i31 *= i32;
      for (k = 0; k + 1 <= i31; k++) {
        i32 = input->size[0];
        slice1_idx_0 = input->size[0];
        input_idx_1 = input->size[1];
        slice1->data[k] = input->data[(k % i32 + slice1_idx_0 * i) +
          slice1_idx_0 * input_idx_1 * (k / i32)];
      }
    }

    e_fillimg(slice1, r12);
    k = output2->size[0];
    slice1_idx_0 = r12->size[0];
    j = r12->size[0];
    loop_ub = j - 1;
    for (i31 = 0; i31 <= loop_ub; i31++) {
      output2->data[i31 + k * i] = r12->data[i31 + slice1_idx_0 * i];
    }
  }

  d_emxInit_real_T(&output3, 7);
  i31 = output3->size[0] * output3->size[1] * output3->size[2] * output3->size[3]
    * output3->size[4] * output3->size[5] * output3->size[6];
  output3->size[0] = input->size[0];
  output3->size[1] = input->size[1];
  output3->size[2] = input->size[2];
  output3->size[3] = input->size[3];
  output3->size[4] = input->size[4];
  output3->size[5] = input->size[5];
  output3->size[6] = input->size[6];
  emxEnsureCapacity((emxArray__common *)output3, i31, (int)sizeof(double));
  loop_ub = input->size[0] * input->size[1] * input->size[2] * input->size[3] *
    input->size[4] * input->size[5] * input->size[6];
  for (i31 = 0; i31 < loop_ub; i31++) {
    output3->data[i31] = input->data[i31];
  }

  i = 0;
  b_emxInit_real_T(&b_slice1, 2);
  while (i <= (int)dim[2] - 1) {
    i31 = input->size[0];
    i32 = input->size[1];
    j = slice1->size[0] * slice1->size[1];
    slice1->size[0] = i31;
    slice1->size[1] = i32;
    emxEnsureCapacity((emxArray__common *)slice1, j, (int)sizeof(double));
    i31 = input->size[0];
    i32 = input->size[1];
    i31 *= i32;
    for (k = 0; k + 1 <= i31; k++) {
      slice1_idx_0 = input->size[0];
      input_idx_1 = input->size[1];
      loop_ub = input->size[0];
      b_loop_ub = input->size[1];
      i32 = r12->size[0] * r12->size[1];
      r12->size[0] = loop_ub;
      r12->size[1] = b_loop_ub;
      emxEnsureCapacity((emxArray__common *)r12, i32, (int)sizeof(double));
      for (i32 = 0; i32 < b_loop_ub; i32++) {
        for (j = 0; j < loop_ub; j++) {
          r12->data[j + r12->size[0] * i32] = input->data[(j + slice1_idx_0 *
            i32) + slice1_idx_0 * input_idx_1 * i];
        }
      }

      i32 = r12->size[0];
      slice1->data[k] = r12->data[k % r12->size[0] + r12->size[0] * (k / i32)];
    }

    i31 = b_slice1->size[0] * b_slice1->size[1];
    b_slice1->size[0] = slice1->size[0];
    b_slice1->size[1] = slice1->size[1];
    emxEnsureCapacity((emxArray__common *)b_slice1, i31, (int)sizeof(double));
    loop_ub = slice1->size[0] * slice1->size[1];
    for (i31 = 0; i31 < loop_ub; i31++) {
      b_slice1->data[i31] = slice1->data[i31];
    }

    e_fillimg(b_slice1, slice1);
    j = output3->size[0];
    k = output3->size[1];
    slice1_idx_0 = slice1->size[0];
    loop_ub = slice1->size[0] - 1;
    b_loop_ub = slice1->size[1] - 1;
    for (i31 = 0; i31 <= b_loop_ub; i31++) {
      for (i32 = 0; i32 <= loop_ub; i32++) {
        output3->data[(i32 + j * i31) + j * k * i] = slice1->data[i32 +
          slice1_idx_0 * i31];
      }
    }

    i++;
  }

  emxFree_real_T(&b_slice1);
  emxFree_real_T(&r12);
  emxFree_real_T(&slice1);
  i31 = output->size[0] * output->size[1] * output->size[2] * output->size[3] *
    output->size[4] * output->size[5] * output->size[6];
  output->size[0] = output1->size[0];
  output->size[1] = output1->size[1];
  output->size[2] = output1->size[2];
  output->size[3] = output1->size[3];
  output->size[4] = output1->size[4];
  output->size[5] = output1->size[5];
  output->size[6] = output1->size[6];
  emxEnsureCapacity((emxArray__common *)output, i31, (int)sizeof(boolean_T));
  loop_ub = output1->size[0] * output1->size[1] * output1->size[2] *
    output1->size[3] * output1->size[4] * output1->size[5] * output1->size[6];
  for (i31 = 0; i31 < loop_ub; i31++) {
    output->data[i31] = ((output1->data[i31] != 0.0) || (output2->data[i31] !=
      0.0) || (output3->data[i31] != 0.0));
  }

  emxFree_real_T(&output3);
  emxFree_real_T(&output2);
  emxFree_real_T(&output1);
}

//
// Arguments    : const emxArray_boolean_T *input
//                emxArray_boolean_T *output
// Return Type  : void
//
void c_fillholes(const emxArray_boolean_T *input, emxArray_boolean_T *output)
{
  int i35;
  int loop_ub;
  unsigned int dim[3];
  int i;
  emxArray_boolean_T *slice1;
  emxArray_boolean_T *r13;
  int b_output;
  boolean_T exitg2;
  int j;
  int sqsz[3];
  int c_output;
  emxArray_boolean_T *output2;
  boolean_T exitg1;
  emxArray_boolean_T *output3;
  emxArray_boolean_T *b_slice1;
  int b_loop_ub;
  i35 = output->size[0] * output->size[1] * output->size[2];
  output->size[0] = input->size[0];
  output->size[1] = input->size[1];
  output->size[2] = input->size[2];
  emxEnsureCapacity((emxArray__common *)output, i35, (int)sizeof(boolean_T));
  loop_ub = input->size[0] * input->size[1] * input->size[2];
  for (i35 = 0; i35 < loop_ub; i35++) {
    output->data[i35] = input->data[i35];
  }

  for (i35 = 0; i35 < 3; i35++) {
    dim[i35] = (unsigned int)input->size[i35];
  }

  i = 0;
  emxInit_boolean_T(&slice1, 2);
  emxInit_boolean_T(&r13, 2);
  while (i <= (int)dim[0] - 1) {
    b_output = 3;
    exitg2 = false;
    while ((!exitg2) && (b_output > 2)) {
      i35 = input->size[2];
      if (i35 == 1) {
        b_output = 2;
      } else {
        exitg2 = true;
      }
    }

    if (b_output <= 2) {
      i35 = input->size[1];
      j = slice1->size[0] * slice1->size[1];
      slice1->size[0] = 1;
      slice1->size[1] = i35;
      emxEnsureCapacity((emxArray__common *)slice1, j, (int)sizeof(boolean_T));
      i35 = input->size[1];
      j = input->size[2];
      i35 *= j;
      for (b_output = 0; b_output + 1 <= i35; b_output++) {
        j = input->size[1];
        slice1->data[b_output] = input->data[(i + input->size[0] * (b_output % j))
          + input->size[0] * input->size[1] * (b_output / j)];
      }
    } else {
      for (i35 = 0; i35 < 3; i35++) {
        sqsz[i35] = 1;
      }

      j = 0;
      i35 = input->size[1];
      if (i35 != 1) {
        i35 = input->size[1];
        sqsz[0] = i35;
        j = 1;
      }

      i35 = input->size[2];
      if (i35 != 1) {
        i35 = input->size[2];
        sqsz[j] = i35;
      }

      i35 = slice1->size[0] * slice1->size[1];
      slice1->size[0] = sqsz[0];
      slice1->size[1] = sqsz[1];
      emxEnsureCapacity((emxArray__common *)slice1, i35, (int)sizeof(boolean_T));
      i35 = input->size[1];
      j = input->size[2];
      i35 *= j;
      for (b_output = 0; b_output + 1 <= i35; b_output++) {
        j = input->size[1];
        slice1->data[b_output] = input->data[(i + input->size[0] * (b_output % j))
          + input->size[0] * input->size[1] * (b_output / j)];
      }
    }

    f_fillimg(slice1, r13);
    b_output = r13->size[0];
    c_output = r13->size[1];
    loop_ub = c_output - 1;
    for (i35 = 0; i35 <= loop_ub; i35++) {
      output->data[i + output->size[0] * i35] = r13->data[i + b_output * i35];
    }

    i++;
  }

  c_emxInit_boolean_T(&output2, 3);
  i35 = output2->size[0] * output2->size[1] * output2->size[2];
  output2->size[0] = input->size[0];
  output2->size[1] = input->size[1];
  output2->size[2] = input->size[2];
  emxEnsureCapacity((emxArray__common *)output2, i35, (int)sizeof(boolean_T));
  loop_ub = input->size[0] * input->size[1] * input->size[2];
  for (i35 = 0; i35 < loop_ub; i35++) {
    output2->data[i35] = input->data[i35];
  }

  for (i = 0; i < (int)dim[1]; i++) {
    b_output = 3;
    exitg1 = false;
    while ((!exitg1) && (b_output > 2)) {
      i35 = input->size[2];
      if (i35 == 1) {
        b_output = 2;
      } else {
        exitg1 = true;
      }
    }

    if (b_output <= 2) {
      i35 = input->size[0];
      j = slice1->size[0] * slice1->size[1];
      slice1->size[0] = i35;
      slice1->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)slice1, j, (int)sizeof(boolean_T));
      i35 = input->size[0];
      j = input->size[2];
      i35 *= j;
      for (b_output = 0; b_output + 1 <= i35; b_output++) {
        j = input->size[0];
        slice1->data[b_output] = input->data[(b_output % j + input->size[0] * i)
          + input->size[0] * input->size[1] * (b_output / j)];
      }
    } else {
      for (i35 = 0; i35 < 3; i35++) {
        sqsz[i35] = 1;
      }

      j = 0;
      i35 = input->size[0];
      if (i35 != 1) {
        i35 = input->size[0];
        sqsz[0] = i35;
        j = 1;
      }

      i35 = input->size[2];
      if (i35 != 1) {
        i35 = input->size[2];
        sqsz[j] = i35;
      }

      i35 = slice1->size[0] * slice1->size[1];
      slice1->size[0] = sqsz[0];
      slice1->size[1] = sqsz[1];
      emxEnsureCapacity((emxArray__common *)slice1, i35, (int)sizeof(boolean_T));
      i35 = input->size[0];
      j = input->size[2];
      i35 *= j;
      for (b_output = 0; b_output + 1 <= i35; b_output++) {
        j = input->size[0];
        slice1->data[b_output] = input->data[(b_output % j + input->size[0] * i)
          + input->size[0] * input->size[1] * (b_output / j)];
      }
    }

    f_fillimg(slice1, r13);
    b_output = r13->size[0];
    c_output = r13->size[0];
    loop_ub = c_output - 1;
    for (i35 = 0; i35 <= loop_ub; i35++) {
      output2->data[i35 + output2->size[0] * i] = r13->data[i35 + b_output * i];
    }
  }

  c_emxInit_boolean_T(&output3, 3);
  i35 = output3->size[0] * output3->size[1] * output3->size[2];
  output3->size[0] = input->size[0];
  output3->size[1] = input->size[1];
  output3->size[2] = input->size[2];
  emxEnsureCapacity((emxArray__common *)output3, i35, (int)sizeof(boolean_T));
  loop_ub = input->size[0] * input->size[1] * input->size[2];
  for (i35 = 0; i35 < loop_ub; i35++) {
    output3->data[i35] = input->data[i35];
  }

  i = 0;
  emxInit_boolean_T(&b_slice1, 2);
  while (i <= (int)dim[2] - 1) {
    i35 = input->size[0];
    j = input->size[1];
    c_output = slice1->size[0] * slice1->size[1];
    slice1->size[0] = i35;
    slice1->size[1] = j;
    emxEnsureCapacity((emxArray__common *)slice1, c_output, (int)sizeof
                      (boolean_T));
    i35 = input->size[0];
    j = input->size[1];
    i35 *= j;
    for (b_output = 0; b_output + 1 <= i35; b_output++) {
      loop_ub = input->size[0];
      b_loop_ub = input->size[1];
      j = r13->size[0] * r13->size[1];
      r13->size[0] = loop_ub;
      r13->size[1] = b_loop_ub;
      emxEnsureCapacity((emxArray__common *)r13, j, (int)sizeof(boolean_T));
      for (j = 0; j < b_loop_ub; j++) {
        for (c_output = 0; c_output < loop_ub; c_output++) {
          r13->data[c_output + r13->size[0] * j] = input->data[(c_output +
            input->size[0] * j) + input->size[0] * input->size[1] * i];
        }
      }

      j = r13->size[0];
      slice1->data[b_output] = r13->data[b_output % r13->size[0] + r13->size[0] *
        (b_output / j)];
    }

    i35 = b_slice1->size[0] * b_slice1->size[1];
    b_slice1->size[0] = slice1->size[0];
    b_slice1->size[1] = slice1->size[1];
    emxEnsureCapacity((emxArray__common *)b_slice1, i35, (int)sizeof(boolean_T));
    loop_ub = slice1->size[0] * slice1->size[1];
    for (i35 = 0; i35 < loop_ub; i35++) {
      b_slice1->data[i35] = slice1->data[i35];
    }

    f_fillimg(b_slice1, slice1);
    c_output = slice1->size[0];
    loop_ub = slice1->size[0] - 1;
    b_loop_ub = slice1->size[1] - 1;
    for (i35 = 0; i35 <= b_loop_ub; i35++) {
      for (j = 0; j <= loop_ub; j++) {
        output3->data[(j + output3->size[0] * i35) + output3->size[0] *
          output3->size[1] * i] = slice1->data[j + c_output * i35];
      }
    }

    i++;
  }

  emxFree_boolean_T(&b_slice1);
  emxFree_boolean_T(&r13);
  emxFree_boolean_T(&slice1);
  i35 = output->size[0] * output->size[1] * output->size[2];
  emxEnsureCapacity((emxArray__common *)output, i35, (int)sizeof(boolean_T));
  j = output->size[0];
  c_output = output->size[1];
  b_output = output->size[2];
  loop_ub = j * c_output * b_output;
  for (i35 = 0; i35 < loop_ub; i35++) {
    output->data[i35] = (output->data[i35] || output2->data[i35] ||
                         output3->data[i35]);
  }

  emxFree_boolean_T(&output3);
  emxFree_boolean_T(&output2);
}

//
// Arguments    : const double input[8530021]
//                boolean_T output[8530021]
// Return Type  : void
//
void fillholes(const double input[8530021], boolean_T output[8530021])
{
  emxArray_real_T *a3;
  emxArray_real_T *b_a3;
  static double output1[8530021];
  int i;
  static double slice1[44197];
  int k;
  int a3_idx_0;
  int i7;
  int i8;
  emxArray_real_T *c_a3;
  static double output2[8530021];
  static double b_slice1[37249];
  emxArray_real_T *d_a3;
  static double output3[8530021];
  int loop_ub;
  b_emxInit_real_T(&a3, 2);
  c_emxInit_real_T(&b_a3, 3);
  for (i = 0; i < 193; i++) {
    for (k = 0; k < 44197; k++) {
      slice1[k] = input[(i + 193 * (k % 229)) + 44197 * (k / 229)];
    }

    b_fillimg(slice1, a3);
    a3_idx_0 = a3->size[0];
    k = a3->size[1];
    i7 = b_a3->size[0] * b_a3->size[1] * b_a3->size[2];
    b_a3->size[0] = 1;
    b_a3->size[1] = k;
    b_a3->size[2] = 1;
    emxEnsureCapacity((emxArray__common *)b_a3, i7, (int)sizeof(double));
    for (i7 = 0; i7 < k; i7++) {
      b_a3->data[b_a3->size[0] * i7] = a3->data[i + a3_idx_0 * i7];
    }

    for (i7 = 0; i7 < 193; i7++) {
      for (i8 = 0; i8 < 229; i8++) {
        output1[(i + 193 * i8) + 44197 * i7] = b_a3->data[i8 + 229 * i7];
      }
    }
  }

  emxFree_real_T(&b_a3);
  c_emxInit_real_T(&c_a3, 3);
  for (i = 0; i < 229; i++) {
    for (k = 0; k < 37249; k++) {
      b_slice1[k] = input[(k % 193 + 193 * i) + 44197 * (k / 193)];
    }

    c_fillimg(b_slice1, a3);
    a3_idx_0 = a3->size[0];
    k = a3->size[0];
    i7 = c_a3->size[0] * c_a3->size[1] * c_a3->size[2];
    c_a3->size[0] = k;
    c_a3->size[1] = 1;
    c_a3->size[2] = 1;
    emxEnsureCapacity((emxArray__common *)c_a3, i7, (int)sizeof(double));
    for (i7 = 0; i7 < k; i7++) {
      c_a3->data[i7] = a3->data[i7 + a3_idx_0 * i];
    }

    for (i7 = 0; i7 < 193; i7++) {
      for (i8 = 0; i8 < 193; i8++) {
        output2[(i8 + 193 * i) + 44197 * i7] = c_a3->data[i8 + 193 * i7];
      }
    }
  }

  emxFree_real_T(&c_a3);
  b_emxInit_real_T(&d_a3, 2);
  for (i = 0; i < 193; i++) {
    for (k = 0; k < 44197; k++) {
      slice1[k] = input[(k % 193 + 193 * (k / 193)) + 44197 * i];
    }

    d_fillimg(slice1, a3);
    a3_idx_0 = a3->size[0];
    k = a3->size[0];
    loop_ub = a3->size[1];
    i7 = d_a3->size[0] * d_a3->size[1];
    d_a3->size[0] = k;
    d_a3->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)d_a3, i7, (int)sizeof(double));
    for (i7 = 0; i7 < loop_ub; i7++) {
      for (i8 = 0; i8 < k; i8++) {
        d_a3->data[i8 + d_a3->size[0] * i7] = a3->data[i8 + a3_idx_0 * i7];
      }
    }

    for (i7 = 0; i7 < 229; i7++) {
      for (i8 = 0; i8 < 193; i8++) {
        output3[(i8 + 193 * i7) + 44197 * i] = d_a3->data[i8 + 193 * i7];
      }
    }
  }

  emxFree_real_T(&d_a3);
  emxFree_real_T(&a3);
  for (i7 = 0; i7 < 8530021; i7++) {
    output[i7] = ((output1[i7] != 0.0) || (output2[i7] != 0.0) || (output3[i7]
      != 0.0));
  }
}
