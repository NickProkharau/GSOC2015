
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
#include "eml_sort.h"
#include "Segmentation1_emxutil.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : const emxArray_real_T *x
//                emxArray_real_T *y
// Return Type  : void
//
void eml_sort(const emxArray_real_T *x, emxArray_real_T *y)
{
  emxArray_real_T *vwork;
  unsigned int unnamed_idx_0;
  int i2;
  int ix;
  int b_i2;
  emxArray_int32_T *iidx;
  emxArray_int32_T *idx0;
  int i;
  int i1;
  int k;
  boolean_T p;
  int j;
  int pEnd;
  int b_p;
  int q;
  int qEnd;
  int kEnd;
  emxInit_real_T(&vwork, 1);
  unnamed_idx_0 = (unsigned int)x->size[0];
  i2 = vwork->size[0];
  vwork->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)vwork, i2, (int)sizeof(double));
  for (i2 = 0; i2 < 2; i2++) {
    ix = y->size[0] * y->size[1];
    y->size[i2] = x->size[i2];
    emxEnsureCapacity((emxArray__common *)y, ix, (int)sizeof(double));
  }

  b_i2 = 1;
  b_emxInit_int32_T(&iidx, 1);
  b_emxInit_int32_T(&idx0, 1);
  for (i = 0; i < 3; i++) {
    i1 = b_i2;
    b_i2 += x->size[0];
    ix = i1;
    for (k = 0; k < x->size[0]; k++) {
      vwork->data[k] = x->data[ix - 1];
      ix++;
    }

    unnamed_idx_0 = (unsigned int)vwork->size[0];
    i2 = iidx->size[0];
    iidx->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)iidx, i2, (int)sizeof(int));
    if (vwork->size[0] == 0) {
    } else {
      for (k = 1; k <= vwork->size[0]; k++) {
        iidx->data[k - 1] = k;
      }

      for (k = 1; k <= vwork->size[0] - 1; k += 2) {
        if ((vwork->data[k - 1] <= vwork->data[k]) || rtIsNaN(vwork->data[k])) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
        } else {
          iidx->data[k - 1] = k + 1;
          iidx->data[k] = k;
        }
      }

      i2 = idx0->size[0];
      idx0->size[0] = vwork->size[0];
      emxEnsureCapacity((emxArray__common *)idx0, i2, (int)sizeof(int));
      ix = vwork->size[0];
      for (i2 = 0; i2 < ix; i2++) {
        idx0->data[i2] = 1;
      }

      ix = 2;
      while (ix < vwork->size[0]) {
        i2 = ix << 1;
        j = 1;
        for (pEnd = 1 + ix; pEnd < vwork->size[0] + 1; pEnd = qEnd + ix) {
          b_p = j;
          q = pEnd - 1;
          qEnd = j + i2;
          if (qEnd > vwork->size[0] + 1) {
            qEnd = vwork->size[0] + 1;
          }

          k = 0;
          kEnd = qEnd - j;
          while (k + 1 <= kEnd) {
            if ((vwork->data[iidx->data[b_p - 1] - 1] <= vwork->data[iidx->
                 data[q] - 1]) || rtIsNaN(vwork->data[iidx->data[q] - 1])) {
              p = true;
            } else {
              p = false;
            }

            if (p) {
              idx0->data[k] = iidx->data[b_p - 1];
              b_p++;
              if (b_p == pEnd) {
                while (q + 1 < qEnd) {
                  k++;
                  idx0->data[k] = iidx->data[q];
                  q++;
                }
              }
            } else {
              idx0->data[k] = iidx->data[q];
              q++;
              if (q + 1 == qEnd) {
                while (b_p < pEnd) {
                  k++;
                  idx0->data[k] = iidx->data[b_p - 1];
                  b_p++;
                }
              }
            }

            k++;
          }

          for (k = 0; k + 1 <= kEnd; k++) {
            iidx->data[(j + k) - 1] = idx0->data[k];
          }

          j = qEnd;
        }

        ix = i2;
      }
    }

    ix = i1;
    for (k = 0; k < x->size[0]; k++) {
      y->data[ix - 1] = vwork->data[iidx->data[k] - 1];
      ix++;
    }
  }

  emxFree_int32_T(&idx0);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}


