
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
#include "imreconstruct.h"
#include "repmat.h"
#include "padarray.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : const double I[44197]
//                emxArray_real_T *I2
// Return Type  : void
//
void b_fillimg(const double I[44197], emxArray_real_T *I2)
{
  unsigned char sizeB[2];
  int k;
  emxArray_real_T *I3;
  emxArray_real_T *mask;
  int b_sizeB[2];
  int i9;
  int i10;
  int j;
  emxArray_real_T *marker;
  for (k = 0; k < 2; k++) {
    sizeB[k] = (unsigned char)(-36 * k + 231);
  }

  b_emxInit_real_T(&I3, 2);
  b_emxInit_real_T(&mask, 2);
  b_sizeB[0] = sizeB[0];
  b_sizeB[1] = sizeB[1];
  repmat(b_sizeB, I3);
  for (i9 = 0; i9 < 2; i9++) {
    i10 = mask->size[0] * mask->size[1];
    mask->size[i9] = I3->size[i9];
    emxEnsureCapacity((emxArray__common *)mask, i10, (int)sizeof(double));
  }

  for (k = 0; k < I3->size[0]; k++) {
    mask->data[k] = rtMinusInf;
  }

  i9 = mask->size[1];
  for (j = 195; j <= i9; j++) {
    i10 = mask->size[0];
    for (k = 0; k < i10; k++) {
      mask->data[k + mask->size[0] * (j - 1)] = rtMinusInf;
    }
  }

  for (j = 0; j < 193; j++) {
    mask->data[mask->size[0] * (j + 1)] = rtMinusInf;
  }

  for (j = 0; j < 193; j++) {
    i9 = mask->size[0];
    for (k = 231; k <= i9; k++) {
      mask->data[(k + mask->size[0] * (j + 1)) - 1] = rtMinusInf;
    }
  }

  for (i9 = 0; i9 < 193; i9++) {
    for (i10 = 0; i10 < 229; i10++) {
      mask->data[(i10 + mask->size[0] * (1 + i9)) + 1] = I[i10 + 229 * i9];
    }
  }

  i9 = mask->size[0] * mask->size[1];
  emxEnsureCapacity((emxArray__common *)mask, i9, (int)sizeof(double));
  k = mask->size[0];
  j = mask->size[1];
  k *= j;
  for (i9 = 0; i9 < k; i9++) {
    mask->data[i9] = 1.0 - mask->data[i9];
  }

  b_emxInit_real_T(&marker, 2);
  i9 = marker->size[0] * marker->size[1];
  marker->size[0] = mask->size[0];
  marker->size[1] = mask->size[1];
  emxEnsureCapacity((emxArray__common *)marker, i9, (int)sizeof(double));
  k = mask->size[0] * mask->size[1];
  for (i9 = 0; i9 < k; i9++) {
    marker->data[i9] = mask->data[i9];
  }
  for (k = 0; k <= mask->size[0] - 3; k++) {
    i9 = marker->size[1];
    for (j = 0; j <= i9 - 3; j++) {
      marker->data[(k + marker->size[0] * (j + 1)) + 1] = rtMinusInf;
    }
  }

  i9 = I2->size[0] * I2->size[1];
  I2->size[0] = marker->size[0];
  I2->size[1] = marker->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i9, (int)sizeof(double));
  k = marker->size[0] * marker->size[1];
  for (i9 = 0; i9 < k; i9++) {
    I2->data[i9] = marker->data[i9];
  }

  imreconstruct(I2, mask);
  i9 = I2->size[0] * I2->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i9, (int)sizeof(double));
  k = I2->size[0];
  j = I2->size[1];
  k *= j;
  emxFree_real_T(&mask);
  for (i9 = 0; i9 < k; i9++) {
    I2->data[i9] = 1.0 - I2->data[i9];
  }

  k = marker->size[0] - 2;
  i9 = I3->size[0] * I3->size[1];
  I3->size[0] = k;
  emxEnsureCapacity((emxArray__common *)I3, i9, (int)sizeof(double));
  k = marker->size[1] - 2;
  i9 = I3->size[0] * I3->size[1];
  I3->size[1] = k;
  emxEnsureCapacity((emxArray__common *)I3, i9, (int)sizeof(double));
  k = (marker->size[0] - 2) * (marker->size[1] - 2);
  for (i9 = 0; i9 < k; i9++) {
    I3->data[i9] = 0.0;
  }

  for (k = 0; k <= marker->size[0] - 3; k++) {
    for (j = 0; j <= marker->size[1] - 3; j++) {
      I3->data[k + I3->size[0] * j] = I2->data[(k + I2->size[0] * (j + 1)) + 1];
    }
  }

  emxFree_real_T(&marker);
  i9 = I2->size[0] * I2->size[1];
  I2->size[0] = I3->size[0];
  I2->size[1] = I3->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i9, (int)sizeof(double));
  k = I3->size[0] * I3->size[1];
  for (i9 = 0; i9 < k; i9++) {
    I2->data[i9] = I3->data[i9];
  }

  emxFree_real_T(&I3);
}

//
// Arguments    : const double I[37249]
//                emxArray_real_T *I2
// Return Type  : void
//
void c_fillimg(const double I[37249], emxArray_real_T *I2)
{
  int iv0[2];
  emxArray_real_T *I3;
  emxArray_real_T *mask;
  int i12;
  int i13;
  int i;
  int j;
  emxArray_real_T *marker;
  iv0[0] = 195;
  iv0[1] = 195;
  b_emxInit_real_T(&I3, 2);
  b_emxInit_real_T(&mask, 2);
  repmat(iv0, I3);
  for (i12 = 0; i12 < 2; i12++) {
    i13 = mask->size[0] * mask->size[1];
    mask->size[i12] = I3->size[i12];
    emxEnsureCapacity((emxArray__common *)mask, i13, (int)sizeof(double));
  }

  for (i = 0; i < I3->size[0]; i++) {
    mask->data[i] = rtMinusInf;
  }

  i12 = mask->size[1];
  for (j = 195; j <= i12; j++) {
    i13 = mask->size[0];
    for (i = 0; i < i13; i++) {
      mask->data[i + mask->size[0] * (j - 1)] = rtMinusInf;
    }
  }

  for (j = 0; j < 193; j++) {
    mask->data[mask->size[0] * (j + 1)] = rtMinusInf;
  }

  for (j = 0; j < 193; j++) {
    i12 = mask->size[0];
    for (i = 195; i <= i12; i++) {
      mask->data[(i + mask->size[0] * (j + 1)) - 1] = rtMinusInf;
    }
  }

  for (i12 = 0; i12 < 193; i12++) {
    for (i13 = 0; i13 < 193; i13++) {
      mask->data[(i13 + mask->size[0] * (1 + i12)) + 1] = I[i13 + 193 * i12];
    }
  }

  i12 = mask->size[0] * mask->size[1];
  emxEnsureCapacity((emxArray__common *)mask, i12, (int)sizeof(double));
  i = mask->size[0];
  j = mask->size[1];
  i *= j;
  for (i12 = 0; i12 < i; i12++) {
    mask->data[i12] = 1.0 - mask->data[i12];
  }

  b_emxInit_real_T(&marker, 2);
  i12 = marker->size[0] * marker->size[1];
  marker->size[0] = mask->size[0];
  marker->size[1] = mask->size[1];
  emxEnsureCapacity((emxArray__common *)marker, i12, (int)sizeof(double));
  i = mask->size[0] * mask->size[1];
  for (i12 = 0; i12 < i; i12++) {
    marker->data[i12] = mask->data[i12];
  }


  for (i = 0; i <= mask->size[0] - 3; i++) {
    i12 = marker->size[1];
    for (j = 0; j <= i12 - 3; j++) {
      marker->data[(i + marker->size[0] * (j + 1)) + 1] = rtMinusInf;
    }
  }

  i12 = I2->size[0] * I2->size[1];
  I2->size[0] = marker->size[0];
  I2->size[1] = marker->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i12, (int)sizeof(double));
  i = marker->size[0] * marker->size[1];
  for (i12 = 0; i12 < i; i12++) {
    I2->data[i12] = marker->data[i12];
  }

  imreconstruct(I2, mask);
  i12 = I2->size[0] * I2->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i12, (int)sizeof(double));
  i = I2->size[0];
  j = I2->size[1];
  i *= j;
  emxFree_real_T(&mask);
  for (i12 = 0; i12 < i; i12++) {
    I2->data[i12] = 1.0 - I2->data[i12];
  }

  i = marker->size[0] - 2;
  i12 = I3->size[0] * I3->size[1];
  I3->size[0] = i;
  emxEnsureCapacity((emxArray__common *)I3, i12, (int)sizeof(double));
  i = marker->size[1] - 2;
  i12 = I3->size[0] * I3->size[1];
  I3->size[1] = i;
  emxEnsureCapacity((emxArray__common *)I3, i12, (int)sizeof(double));
  i = (marker->size[0] - 2) * (marker->size[1] - 2);
  for (i12 = 0; i12 < i; i12++) {
    I3->data[i12] = 0.0;
  }

  for (i = 0; i <= marker->size[0] - 3; i++) {
    for (j = 0; j <= marker->size[1] - 3; j++) {
      I3->data[i + I3->size[0] * j] = I2->data[(i + I2->size[0] * (j + 1)) + 1];
    }
  }

  emxFree_real_T(&marker);
  i12 = I2->size[0] * I2->size[1];
  I2->size[0] = I3->size[0];
  I2->size[1] = I3->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i12, (int)sizeof(double));
  i = I3->size[0] * I3->size[1];
  for (i12 = 0; i12 < i; i12++) {
    I2->data[i12] = I3->data[i12];
  }

  emxFree_real_T(&I3);
}

//
// Arguments    : const double I[44197]
//                emxArray_real_T *I2
// Return Type  : void
//
void d_fillimg(const double I[44197], emxArray_real_T *I2)
{
  unsigned char sizeB[2];
  int k;
  emxArray_real_T *I3;
  emxArray_real_T *mask;
  int b_sizeB[2];
  int i14;
  int i15;
  int j;
  emxArray_real_T *marker;
  for (k = 0; k < 2; k++) {
    sizeB[k] = (unsigned char)(36 * k + 195);
  }

  b_emxInit_real_T(&I3, 2);
  b_emxInit_real_T(&mask, 2);
  b_sizeB[0] = sizeB[0];
  b_sizeB[1] = sizeB[1];
  repmat(b_sizeB, I3);
  for (i14 = 0; i14 < 2; i14++) {
    i15 = mask->size[0] * mask->size[1];
    mask->size[i14] = I3->size[i14];
    emxEnsureCapacity((emxArray__common *)mask, i15, (int)sizeof(double));
  }

  for (k = 0; k < I3->size[0]; k++) {
    mask->data[k] = rtMinusInf;
  }

  i14 = mask->size[1];
  for (j = 231; j <= i14; j++) {
    i15 = mask->size[0];
    for (k = 0; k < i15; k++) {
      mask->data[k + mask->size[0] * (j - 1)] = rtMinusInf;
    }
  }

  for (j = 0; j < 229; j++) {
    mask->data[mask->size[0] * (j + 1)] = rtMinusInf;
  }

  for (j = 0; j < 229; j++) {
    i14 = mask->size[0];
    for (k = 195; k <= i14; k++) {
      mask->data[(k + mask->size[0] * (j + 1)) - 1] = rtMinusInf;
    }
  }

  for (i14 = 0; i14 < 229; i14++) {
    for (i15 = 0; i15 < 193; i15++) {
      mask->data[(i15 + mask->size[0] * (1 + i14)) + 1] = I[i15 + 193 * i14];
    }
  }

  i14 = mask->size[0] * mask->size[1];
  emxEnsureCapacity((emxArray__common *)mask, i14, (int)sizeof(double));
  k = mask->size[0];
  j = mask->size[1];
  k *= j;
  for (i14 = 0; i14 < k; i14++) {
    mask->data[i14] = 1.0 - mask->data[i14];
  }

  b_emxInit_real_T(&marker, 2);
  i14 = marker->size[0] * marker->size[1];
  marker->size[0] = mask->size[0];
  marker->size[1] = mask->size[1];
  emxEnsureCapacity((emxArray__common *)marker, i14, (int)sizeof(double));
  k = mask->size[0] * mask->size[1];
  for (i14 = 0; i14 < k; i14++) {
    marker->data[i14] = mask->data[i14];
  }

  
  for (k = 0; k <= mask->size[0] - 3; k++) {
    i14 = marker->size[1];
    for (j = 0; j <= i14 - 3; j++) {
      marker->data[(k + marker->size[0] * (j + 1)) + 1] = rtMinusInf;
    }
  }

  i14 = I2->size[0] * I2->size[1];
  I2->size[0] = marker->size[0];
  I2->size[1] = marker->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i14, (int)sizeof(double));
  k = marker->size[0] * marker->size[1];
  for (i14 = 0; i14 < k; i14++) {
    I2->data[i14] = marker->data[i14];
  }

  imreconstruct(I2, mask);
  i14 = I2->size[0] * I2->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i14, (int)sizeof(double));
  k = I2->size[0];
  j = I2->size[1];
  k *= j;
  emxFree_real_T(&mask);
  for (i14 = 0; i14 < k; i14++) {
    I2->data[i14] = 1.0 - I2->data[i14];
  }

  k = marker->size[0] - 2;
  i14 = I3->size[0] * I3->size[1];
  I3->size[0] = k;
  emxEnsureCapacity((emxArray__common *)I3, i14, (int)sizeof(double));
  k = marker->size[1] - 2;
  i14 = I3->size[0] * I3->size[1];
  I3->size[1] = k;
  emxEnsureCapacity((emxArray__common *)I3, i14, (int)sizeof(double));
  k = (marker->size[0] - 2) * (marker->size[1] - 2);
  for (i14 = 0; i14 < k; i14++) {
    I3->data[i14] = 0.0;
  }

  for (k = 0; k <= marker->size[0] - 3; k++) {
    for (j = 0; j <= marker->size[1] - 3; j++) {
      I3->data[k + I3->size[0] * j] = I2->data[(k + I2->size[0] * (j + 1)) + 1];
    }
  }

  emxFree_real_T(&marker);
  i14 = I2->size[0] * I2->size[1];
  I2->size[0] = I3->size[0];
  I2->size[1] = I3->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i14, (int)sizeof(double));
  k = I3->size[0] * I3->size[1];
  for (i14 = 0; i14 < k; i14++) {
    I2->data[i14] = I3->data[i14];
  }

  emxFree_real_T(&I3);
}

//
// Arguments    : const emxArray_real_T *I
//                emxArray_real_T *I2
// Return Type  : void
//
void e_fillimg(const emxArray_real_T *I, emxArray_real_T *I2)
{
  emxArray_real_T *mask;
  int i33;
  int b_mask;
  int c_mask;
  emxArray_real_T *marker;
  b_emxInit_real_T(&mask, 2);
  padarray(I, mask);
  i33 = mask->size[0] * mask->size[1];
  emxEnsureCapacity((emxArray__common *)mask, i33, (int)sizeof(double));
  b_mask = mask->size[0];
  c_mask = mask->size[1];
  b_mask *= c_mask;
  for (i33 = 0; i33 < b_mask; i33++) {
    mask->data[i33] = 1.0 - mask->data[i33];
  }

  b_emxInit_real_T(&marker, 2);
  i33 = marker->size[0] * marker->size[1];
  marker->size[0] = mask->size[0];
  marker->size[1] = mask->size[1];
  emxEnsureCapacity((emxArray__common *)marker, i33, (int)sizeof(double));
  b_mask = mask->size[0] * mask->size[1];
  for (i33 = 0; i33 < b_mask; i33++) {
    marker->data[i33] = mask->data[i33];
  }

 
  for (b_mask = 0; b_mask <= mask->size[0] - 3; b_mask++) {
    i33 = marker->size[1];
    for (c_mask = 0; c_mask <= i33 - 3; c_mask++) {
      marker->data[(b_mask + marker->size[0] * (c_mask + 1)) + 1] = rtMinusInf;
    }
  }

  i33 = I2->size[0] * I2->size[1];
  I2->size[0] = marker->size[0];
  I2->size[1] = marker->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i33, (int)sizeof(double));
  b_mask = marker->size[0] * marker->size[1];
  for (i33 = 0; i33 < b_mask; i33++) {
    I2->data[i33] = marker->data[i33];
  }

  imreconstruct(I2, mask);
  i33 = I2->size[0] * I2->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i33, (int)sizeof(double));
  b_mask = I2->size[0];
  c_mask = I2->size[1];
  b_mask *= c_mask;
  for (i33 = 0; i33 < b_mask; i33++) {
    I2->data[i33] = 1.0 - I2->data[i33];
  }

  b_mask = marker->size[0] - 2;
  i33 = mask->size[0] * mask->size[1];
  mask->size[0] = b_mask;
  emxEnsureCapacity((emxArray__common *)mask, i33, (int)sizeof(double));
  b_mask = marker->size[1] - 2;
  i33 = mask->size[0] * mask->size[1];
  mask->size[1] = b_mask;
  emxEnsureCapacity((emxArray__common *)mask, i33, (int)sizeof(double));
  b_mask = (marker->size[0] - 2) * (marker->size[1] - 2);
  for (i33 = 0; i33 < b_mask; i33++) {
    mask->data[i33] = 0.0;
  }

  for (b_mask = 0; b_mask <= marker->size[0] - 3; b_mask++) {
    for (c_mask = 0; c_mask <= marker->size[1] - 3; c_mask++) {
      mask->data[b_mask + mask->size[0] * c_mask] = I2->data[(b_mask + I2->size
        [0] * (c_mask + 1)) + 1];
    }
  }

  emxFree_real_T(&marker);
  i33 = I2->size[0] * I2->size[1];
  I2->size[0] = mask->size[0];
  I2->size[1] = mask->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i33, (int)sizeof(double));
  b_mask = mask->size[0] * mask->size[1];
  for (i33 = 0; i33 < b_mask; i33++) {
    I2->data[i33] = mask->data[i33];
  }

  emxFree_real_T(&mask);
}

//
// Arguments    : const emxArray_boolean_T *I
//                emxArray_boolean_T *I2
// Return Type  : void
//
void f_fillimg(const emxArray_boolean_T *I, emxArray_boolean_T *I2)
{
  emxArray_uint8_T *mask;
  int i36;
  int loop_ub;
  double sizeB[2];
  emxArray_uint8_T *b_mask;
  int c_mask;
  emxArray_uint8_T *marker;
  emxArray_uint8_T *b_I2;
  emxArray_real_T *I3;
  b_emxInit_uint8_T(&mask, 2);
  i36 = mask->size[0] * mask->size[1];
  mask->size[0] = I->size[0];
  mask->size[1] = I->size[1];
  emxEnsureCapacity((emxArray__common *)mask, i36, (int)sizeof(unsigned char));
  loop_ub = I->size[0] * I->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    mask->data[i36] = I->data[i36];
  }

  if ((mask->size[0] == 0) || (mask->size[1] == 0)) {
    for (i36 = 0; i36 < 2; i36++) {
      sizeB[i36] = (double)mask->size[i36] + 2.0;
    }

    i36 = mask->size[0] * mask->size[1];
    mask->size[0] = (int)sizeB[0];
    emxEnsureCapacity((emxArray__common *)mask, i36, (int)sizeof(unsigned char));
    i36 = mask->size[0] * mask->size[1];
    mask->size[1] = (int)sizeB[1];
    emxEnsureCapacity((emxArray__common *)mask, i36, (int)sizeof(unsigned char));
    loop_ub = (int)sizeB[0] * (int)sizeB[1];
    for (i36 = 0; i36 < loop_ub; i36++) {
      mask->data[i36] = 0;
    }
  } else {
    b_emxInit_uint8_T(&b_mask, 2);
    i36 = b_mask->size[0] * b_mask->size[1];
    b_mask->size[0] = mask->size[0];
    b_mask->size[1] = mask->size[1];
    emxEnsureCapacity((emxArray__common *)b_mask, i36, (int)sizeof(unsigned char));
    loop_ub = mask->size[0] * mask->size[1];
    for (i36 = 0; i36 < loop_ub; i36++) {
      b_mask->data[i36] = mask->data[i36];
    }

    b_ConstantPad(b_mask, 7, mask);
    emxFree_uint8_T(&b_mask);
  }

  i36 = mask->size[0] * mask->size[1];
  emxEnsureCapacity((emxArray__common *)mask, i36, (int)sizeof(unsigned char));
  loop_ub = mask->size[0];
  c_mask = mask->size[1];
  loop_ub *= c_mask;
  for (i36 = 0; i36 < loop_ub; i36++) {
    mask->data[i36] = (unsigned char)(255U - mask->data[i36]);
  }

  b_emxInit_uint8_T(&marker, 2);
  i36 = marker->size[0] * marker->size[1];
  marker->size[0] = mask->size[0];
  marker->size[1] = mask->size[1];
  emxEnsureCapacity((emxArray__common *)marker, i36, (int)sizeof(unsigned char));
  loop_ub = mask->size[0] * mask->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    marker->data[i36] = mask->data[i36];
  }


  for (loop_ub = 0; loop_ub <= mask->size[0] - 3; loop_ub++) {
    i36 = marker->size[1];
    for (c_mask = 0; c_mask <= i36 - 3; c_mask++) {
      marker->data[(loop_ub + marker->size[0] * (c_mask + 1)) + 1] = 0;
    }
  }

  for (i36 = 0; i36 < 2; i36++) {
    sizeB[i36] = marker->size[i36];
  }

  b_emxInit_uint8_T(&b_I2, 2);
  i36 = b_I2->size[0] * b_I2->size[1];
  b_I2->size[0] = marker->size[0];
  b_I2->size[1] = marker->size[1];
  emxEnsureCapacity((emxArray__common *)b_I2, i36, (int)sizeof(unsigned char));
  loop_ub = marker->size[0] * marker->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    b_I2->data[i36] = marker->data[i36];
  }

  //ippreconstruct_uint8(&b_I2->data[0], &mask->data[0], sizeB, 2.0, 2.0);
  i36 = b_I2->size[0] * b_I2->size[1];
  emxEnsureCapacity((emxArray__common *)b_I2, i36, (int)sizeof(unsigned char));
  loop_ub = b_I2->size[0];
  c_mask = b_I2->size[1];
  loop_ub *= c_mask;
  emxFree_uint8_T(&mask);
  for (i36 = 0; i36 < loop_ub; i36++) {
    b_I2->data[i36] = (unsigned char)(255U - b_I2->data[i36]);
  }

  b_emxInit_real_T(&I3, 2);
  loop_ub = marker->size[0] - 2;
  i36 = I3->size[0] * I3->size[1];
  I3->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)I3, i36, (int)sizeof(double));
  loop_ub = marker->size[1] - 2;
  i36 = I3->size[0] * I3->size[1];
  I3->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)I3, i36, (int)sizeof(double));
  loop_ub = (marker->size[0] - 2) * (marker->size[1] - 2);
  for (i36 = 0; i36 < loop_ub; i36++) {
    I3->data[i36] = 0.0;
  }

  for (loop_ub = 0; loop_ub <= marker->size[0] - 3; loop_ub++) {
    for (c_mask = 0; c_mask <= marker->size[1] - 3; c_mask++) {
      I3->data[loop_ub + I3->size[0] * c_mask] = b_I2->data[(loop_ub +
        b_I2->size[0] * (c_mask + 1)) + 1];
    }
  }

  emxFree_uint8_T(&b_I2);
  emxFree_uint8_T(&marker);
  i36 = I2->size[0] * I2->size[1];
  I2->size[0] = I3->size[0];
  I2->size[1] = I3->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i36, (int)sizeof(boolean_T));
  loop_ub = I3->size[0] * I3->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    I2->data[i36] = (I3->data[i36] != 0.0);
  }

  emxFree_real_T(&I3);
}

//
// Arguments    : const emxArray_real_T *I
//                emxArray_real_T *I2
// Return Type  : void
//
void fillimg(const emxArray_real_T *I, emxArray_real_T *I2)
{
  emxArray_real_T *mask;
  int i16;
  int b_mask;
  int c_mask;
  emxArray_real_T *marker;
  b_emxInit_real_T(&mask, 2);
  padarray(I, mask);
  i16 = mask->size[0] * mask->size[1];
  emxEnsureCapacity((emxArray__common *)mask, i16, (int)sizeof(double));
  b_mask = mask->size[0];
  c_mask = mask->size[1];
  b_mask *= c_mask;
  for (i16 = 0; i16 < b_mask; i16++) {
    mask->data[i16] = 1.0 - mask->data[i16];
  }

  b_emxInit_real_T(&marker, 2);
  i16 = marker->size[0] * marker->size[1];
  marker->size[0] = mask->size[0];
  marker->size[1] = mask->size[1];
  emxEnsureCapacity((emxArray__common *)marker, i16, (int)sizeof(double));
  b_mask = mask->size[0] * mask->size[1];
  for (i16 = 0; i16 < b_mask; i16++) {
    marker->data[i16] = mask->data[i16];
  }

  for (b_mask = 0; b_mask <= mask->size[0] - 3; b_mask++) {
    i16 = marker->size[1];
    for (c_mask = 0; c_mask <= i16 - 3; c_mask++) {
      marker->data[(b_mask + marker->size[0] * (c_mask + 1)) + 1] = rtMinusInf;
    }
  }

  i16 = I2->size[0] * I2->size[1];
  I2->size[0] = marker->size[0];
  I2->size[1] = marker->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i16, (int)sizeof(double));
  b_mask = marker->size[0] * marker->size[1];
  for (i16 = 0; i16 < b_mask; i16++) {
    I2->data[i16] = marker->data[i16];
  }

  imreconstruct(I2, mask);
  i16 = I2->size[0] * I2->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i16, (int)sizeof(double));
  b_mask = I2->size[0];
  c_mask = I2->size[1];
  b_mask *= c_mask;
  for (i16 = 0; i16 < b_mask; i16++) {
    I2->data[i16] = 1.0 - I2->data[i16];
  }

  b_mask = marker->size[0] - 2;
  i16 = mask->size[0] * mask->size[1];
  mask->size[0] = b_mask;
  emxEnsureCapacity((emxArray__common *)mask, i16, (int)sizeof(double));
  b_mask = marker->size[1] - 2;
  i16 = mask->size[0] * mask->size[1];
  mask->size[1] = b_mask;
  emxEnsureCapacity((emxArray__common *)mask, i16, (int)sizeof(double));
  b_mask = (marker->size[0] - 2) * (marker->size[1] - 2);
  for (i16 = 0; i16 < b_mask; i16++) {
    mask->data[i16] = 0.0;
  }

  for (b_mask = 0; b_mask <= marker->size[0] - 3; b_mask++) {
    for (c_mask = 0; c_mask <= marker->size[1] - 3; c_mask++) {
      mask->data[b_mask + mask->size[0] * c_mask] = I2->data[(b_mask + I2->size
        [0] * (c_mask + 1)) + 1];
    }
  }

  emxFree_real_T(&marker);
  i16 = I2->size[0] * I2->size[1];
  I2->size[0] = mask->size[0];
  I2->size[1] = mask->size[1];
  emxEnsureCapacity((emxArray__common *)I2, i16, (int)sizeof(double));
  b_mask = mask->size[0] * mask->size[1];
  for (i16 = 0; i16 < b_mask; i16++) {
    I2->data[i16] = mask->data[i16];
  }

  emxFree_real_T(&mask);
}
