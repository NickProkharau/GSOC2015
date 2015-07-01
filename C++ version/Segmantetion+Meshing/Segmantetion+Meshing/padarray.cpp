
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
#include "padarray.h"
#include "Segmentation1_emxutil.h"
#include "repmat.h"
#include <stdio.h>


// Function Declarations
static void ConstantPad(const emxArray_real_T *a, signed char direction,
  emxArray_real_T *b);

// Function Definitions

//
// Arguments    : const emxArray_real_T *a
//                signed char direction
//                emxArray_real_T *b
// Return Type  : void
//
static void ConstantPad(const emxArray_real_T *a, signed char direction,
  emxArray_real_T *b)
{
  int sizeB[2];
  int cdiff;
  int i18;
  emxArray_real_T *r0;
  int b_sizeB[2];
  int absb;
  int ndbl;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  int apnd;
  for (cdiff = 0; cdiff < 2; cdiff++) {
    if (direction == 5) {
      i18 = a->size[cdiff] + 1;
      sizeB[cdiff] = i18;
    } else if (direction == 6) {
      i18 = a->size[cdiff] + 1;
      sizeB[cdiff] = i18;
    } else {
      i18 = a->size[cdiff];
      sizeB[cdiff] = (int)((double)i18 + 2.0);
    }
  }

  b_emxInit_real_T(&r0, 2);
  b_sizeB[0] = sizeB[0];
  b_sizeB[1] = sizeB[1];
  repmat(b_sizeB, r0);
  for (i18 = 0; i18 < 2; i18++) {
    absb = b->size[0] * b->size[1];
    b->size[i18] = r0->size[i18];
    emxEnsureCapacity((emxArray__common *)b, absb, (int)sizeof(double));
  }

  if (direction == 5) {
    for (cdiff = 0; cdiff < r0->size[0]; cdiff++) {
      b->data[cdiff] = rtMinusInf;
    }

    i18 = b->size[1];
    for (ndbl = 2; ndbl <= i18; ndbl++) {
      b->data[b->size[0] * (ndbl - 1)] = rtMinusInf;
    }
  } else if (direction == 7) {
    for (cdiff = 0; cdiff < r0->size[0]; cdiff++) {
      b->data[cdiff] = rtMinusInf;
    }

    i18 = b->size[1];
    for (ndbl = a->size[1] + 1; ndbl + 1 <= i18; ndbl++) {
      absb = b->size[0];
      for (cdiff = 0; cdiff < absb; cdiff++) {
        b->data[cdiff + b->size[0] * ndbl] = rtMinusInf;
      }
    }

    for (ndbl = 2; ndbl <= a->size[1] + 1; ndbl++) {
      b->data[b->size[0] * (ndbl - 1)] = rtMinusInf;
    }

    for (ndbl = 2; ndbl <= a->size[1] + 1; ndbl++) {
      i18 = b->size[0];
      for (cdiff = a->size[0] + 1; cdiff + 1 <= i18; cdiff++) {
        b->data[cdiff + b->size[0] * (ndbl - 1)] = rtMinusInf;
      }
    }
  } else {
    for (ndbl = a->size[1]; ndbl + 1 <= r0->size[1]; ndbl++) {
      i18 = b->size[0];
      for (cdiff = 0; cdiff < i18; cdiff++) {
        b->data[cdiff + b->size[0] * ndbl] = rtMinusInf;
      }
    }

    for (ndbl = 0; ndbl < a->size[1]; ndbl++) {
      i18 = b->size[0];
      for (cdiff = a->size[0]; cdiff + 1 <= i18; cdiff++) {
        b->data[cdiff + b->size[0] * ndbl] = rtMinusInf;
      }
    }
  }

  emxFree_real_T(&r0);
  b_emxInit_real_T(&y, 2);
  b_emxInit_real_T(&b_y, 2);
  if ((direction == 5) || (direction == 7)) {
    if (a->size[0] < 1) {
      absb = -1;
      apnd = 0;
    } else {
      ndbl = (int)floor(((double)a->size[0] - 1.0) + 0.5);
      apnd = ndbl + 1;
      cdiff = (ndbl - a->size[0]) + 1;
      absb = a->size[0];
      if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
        ndbl++;
        apnd = a->size[0];
      } else if (cdiff > 0) {
        apnd = ndbl;
      } else {
        ndbl++;
      }

      absb = ndbl - 1;
    }

    i18 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = absb + 1;
    emxEnsureCapacity((emxArray__common *)y, i18, (int)sizeof(double));
    if (absb + 1 > 0) {
      y->data[0] = 1.0;
      if (absb + 1 > 1) {
        y->data[absb] = apnd;
        ndbl = absb / 2;
        for (cdiff = 1; cdiff < ndbl; cdiff++) {
          y->data[cdiff] = 1.0 + (double)cdiff;
          y->data[absb - cdiff] = apnd - cdiff;
        }

        if (ndbl << 1 == absb) {
          y->data[ndbl] = (1.0 + (double)apnd) / 2.0;
        } else {
          y->data[ndbl] = 1.0 + (double)ndbl;
          y->data[ndbl + 1] = apnd - ndbl;
        }
      }
    }

    if (a->size[1] < 1) {
      absb = -1;
      apnd = 0;
    } else {
      ndbl = (int)floor(((double)a->size[1] - 1.0) + 0.5);
      apnd = ndbl + 1;
      cdiff = (ndbl - a->size[1]) + 1;
      absb = a->size[1];
      if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
        ndbl++;
        apnd = a->size[1];
      } else if (cdiff > 0) {
        apnd = ndbl;
      } else {
        ndbl++;
      }

      absb = ndbl - 1;
    }

    i18 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = absb + 1;
    emxEnsureCapacity((emxArray__common *)b_y, i18, (int)sizeof(double));
    if (absb + 1 > 0) {
      b_y->data[0] = 1.0;
      if (absb + 1 > 1) {
        b_y->data[absb] = apnd;
        ndbl = absb / 2;
        for (cdiff = 1; cdiff < ndbl; cdiff++) {
          b_y->data[cdiff] = 1.0 + (double)cdiff;
          b_y->data[absb - cdiff] = apnd - cdiff;
        }

        if (ndbl << 1 == absb) {
          b_y->data[ndbl] = (1.0 + (double)apnd) / 2.0;
        } else {
          b_y->data[ndbl] = 1.0 + (double)ndbl;
          b_y->data[ndbl + 1] = apnd - ndbl;
        }
      }
    }

    ndbl = a->size[1];
    for (i18 = 0; i18 < ndbl; i18++) {
      cdiff = a->size[0];
      for (absb = 0; absb < cdiff; absb++) {
        b->data[(int)y->data[y->size[0] * absb] + b->size[0] * (int)b_y->
          data[b_y->size[0] * i18]] = a->data[absb + a->size[0] * i18];
      }
    }
  } else {
    ndbl = a->size[1];
    for (i18 = 0; i18 < ndbl; i18++) {
      cdiff = a->size[0];
      for (absb = 0; absb < cdiff; absb++) {
        b->data[absb + b->size[0] * i18] = a->data[absb + a->size[0] * i18];
      }
    }
  }

  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
}

//
// Arguments    : const emxArray_uint8_T *a
//                signed char direction
//                emxArray_uint8_T *b
// Return Type  : void
//
void b_ConstantPad(const emxArray_uint8_T *a, signed char direction,
                   emxArray_uint8_T *b)
{
  int sizeB[2];
  int cdiff;
  int i37;
  emxArray_uint8_T *r14;
  int b_sizeB[2];
  int absb;
  int ndbl;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  int apnd;
  for (cdiff = 0; cdiff < 2; cdiff++) {
    if (direction == 5) {
      i37 = a->size[cdiff] + 1;
      sizeB[cdiff] = i37;
    } else if (direction == 6) {
      i37 = a->size[cdiff] + 1;
      sizeB[cdiff] = i37;
    } else {
      i37 = a->size[cdiff];
      sizeB[cdiff] = (int)((double)i37 + 2.0);
    }
  }

  b_emxInit_uint8_T(&r14, 2);
  b_sizeB[0] = sizeB[0];
  b_sizeB[1] = sizeB[1];
  b_repmat(b_sizeB, r14);
  for (i37 = 0; i37 < 2; i37++) {
    absb = b->size[0] * b->size[1];
    b->size[i37] = r14->size[i37];
    emxEnsureCapacity((emxArray__common *)b, absb, (int)sizeof(unsigned char));
  }

  if (direction == 5) {
    for (cdiff = 0; cdiff < r14->size[0]; cdiff++) {
      b->data[cdiff] = 0;
    }

    i37 = b->size[1];
    for (ndbl = 2; ndbl <= i37; ndbl++) {
      b->data[b->size[0] * (ndbl - 1)] = 0;
    }
  } else if (direction == 7) {
    for (cdiff = 0; cdiff < r14->size[0]; cdiff++) {
      b->data[cdiff] = 0;
    }

    i37 = b->size[1];
    for (ndbl = a->size[1] + 1; ndbl + 1 <= i37; ndbl++) {
      absb = b->size[0];
      for (cdiff = 0; cdiff < absb; cdiff++) {
        b->data[cdiff + b->size[0] * ndbl] = 0;
      }
    }

    for (ndbl = 2; ndbl <= a->size[1] + 1; ndbl++) {
      b->data[b->size[0] * (ndbl - 1)] = 0;
    }

    for (ndbl = 2; ndbl <= a->size[1] + 1; ndbl++) {
      i37 = b->size[0];
      for (cdiff = a->size[0] + 1; cdiff + 1 <= i37; cdiff++) {
        b->data[cdiff + b->size[0] * (ndbl - 1)] = 0;
      }
    }
  } else {
    for (ndbl = a->size[1]; ndbl + 1 <= r14->size[1]; ndbl++) {
      i37 = b->size[0];
      for (cdiff = 0; cdiff < i37; cdiff++) {
        b->data[cdiff + b->size[0] * ndbl] = 0;
      }
    }

    for (ndbl = 0; ndbl < a->size[1]; ndbl++) {
      i37 = b->size[0];
      for (cdiff = a->size[0]; cdiff + 1 <= i37; cdiff++) {
        b->data[cdiff + b->size[0] * ndbl] = 0;
      }
    }
  }

  emxFree_uint8_T(&r14);
  b_emxInit_real_T(&y, 2);
  b_emxInit_real_T(&b_y, 2);
  if ((direction == 5) || (direction == 7)) {
    if (a->size[0] < 1) {
      absb = -1;
      apnd = 0;
    } else {
      ndbl = (int)floor(((double)a->size[0] - 1.0) + 0.5);
      apnd = ndbl + 1;
      cdiff = (ndbl - a->size[0]) + 1;
      absb = a->size[0];
      if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
        ndbl++;
        apnd = a->size[0];
      } else if (cdiff > 0) {
        apnd = ndbl;
      } else {
        ndbl++;
      }

      absb = ndbl - 1;
    }

    i37 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = absb + 1;
    emxEnsureCapacity((emxArray__common *)y, i37, (int)sizeof(double));
    if (absb + 1 > 0) {
      y->data[0] = 1.0;
      if (absb + 1 > 1) {
        y->data[absb] = apnd;
        ndbl = absb / 2;
        for (cdiff = 1; cdiff < ndbl; cdiff++) {
          y->data[cdiff] = 1.0 + (double)cdiff;
          y->data[absb - cdiff] = apnd - cdiff;
        }

        if (ndbl << 1 == absb) {
          y->data[ndbl] = (1.0 + (double)apnd) / 2.0;
        } else {
          y->data[ndbl] = 1.0 + (double)ndbl;
          y->data[ndbl + 1] = apnd - ndbl;
        }
      }
    }

    if (a->size[1] < 1) {
      absb = -1;
      apnd = 0;
    } else {
      ndbl = (int)floor(((double)a->size[1] - 1.0) + 0.5);
      apnd = ndbl + 1;
      cdiff = (ndbl - a->size[1]) + 1;
      absb = a->size[1];
      if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
        ndbl++;
        apnd = a->size[1];
      } else if (cdiff > 0) {
        apnd = ndbl;
      } else {
        ndbl++;
      }

      absb = ndbl - 1;
    }

    i37 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = absb + 1;
    emxEnsureCapacity((emxArray__common *)b_y, i37, (int)sizeof(double));
    if (absb + 1 > 0) {
      b_y->data[0] = 1.0;
      if (absb + 1 > 1) {
        b_y->data[absb] = apnd;
        ndbl = absb / 2;
        for (cdiff = 1; cdiff < ndbl; cdiff++) {
          b_y->data[cdiff] = 1.0 + (double)cdiff;
          b_y->data[absb - cdiff] = apnd - cdiff;
        }

        if (ndbl << 1 == absb) {
          b_y->data[ndbl] = (1.0 + (double)apnd) / 2.0;
        } else {
          b_y->data[ndbl] = 1.0 + (double)ndbl;
          b_y->data[ndbl + 1] = apnd - ndbl;
        }
      }
    }

    ndbl = a->size[1];
    for (i37 = 0; i37 < ndbl; i37++) {
      cdiff = a->size[0];
      for (absb = 0; absb < cdiff; absb++) {
        b->data[(int)y->data[y->size[0] * absb] + b->size[0] * (int)b_y->
          data[b_y->size[0] * i37]] = a->data[absb + a->size[0] * i37];
      }
    }
  } else {
    ndbl = a->size[1];
    for (i37 = 0; i37 < ndbl; i37++) {
      cdiff = a->size[0];
      for (absb = 0; absb < cdiff; absb++) {
        b->data[absb + b->size[0] * i37] = a->data[absb + a->size[0] * i37];
      }
    }
  }

  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
}

//
// Arguments    : const emxArray_real_T *varargin_1
//                emxArray_real_T *b
// Return Type  : void
//
void padarray(const emxArray_real_T *varargin_1, emxArray_real_T *b)
{
  unsigned int sizeB[2];
  int i17;
  int loop_ub;
  if ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0)) {
    for (i17 = 0; i17 < 2; i17++) {
      sizeB[i17] = varargin_1->size[i17] + 2U;
    }

    i17 = b->size[0] * b->size[1];
    b->size[0] = (int)sizeB[0];
    emxEnsureCapacity((emxArray__common *)b, i17, (int)sizeof(double));
    i17 = b->size[0] * b->size[1];
    b->size[1] = (int)sizeB[1];
    emxEnsureCapacity((emxArray__common *)b, i17, (int)sizeof(double));
    loop_ub = (int)sizeB[0] * (int)sizeB[1];
    for (i17 = 0; i17 < loop_ub; i17++) {
      b->data[i17] = rtMinusInf;
    }
  } else {
    ConstantPad(varargin_1, 7, b);
  }
}

