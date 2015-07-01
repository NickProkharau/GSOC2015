
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
#include "norm.h"
#include "Segmentation1_emxutil.h"
#include <stdio.h>


// Function Definitions

//
// sphere constructs a 3D sphere with the specified radius
//  that can be used as structural element in 3D image processing
// Arguments    : double r
//                emxArray_real_T *se
// Return Type  : void
//
void sphere(double r, emxArray_real_T *se)
{
  double dim_idx_0;
  double dim_idx_1;
  double dim_idx_2;
  int i;
  int loop_ub;
  int k;
  double dv5[3];
  dim_idx_0 = 2.0 * r + 1.0;
  dim_idx_1 = 2.0 * r + 1.0;
  dim_idx_2 = 2.0 * r + 1.0;
  i = se->size[0] * se->size[1] * se->size[2];
  se->size[0] = (int)dim_idx_0;
  emxEnsureCapacity((emxArray__common *)se, i, (int)sizeof(double));
  i = se->size[0] * se->size[1] * se->size[2];
  se->size[1] = (int)dim_idx_1;
  emxEnsureCapacity((emxArray__common *)se, i, (int)sizeof(double));
  i = se->size[0] * se->size[1] * se->size[2];
  se->size[2] = (int)dim_idx_2;
  emxEnsureCapacity((emxArray__common *)se, i, (int)sizeof(double));
  loop_ub = (int)dim_idx_0 * (int)dim_idx_1 * (int)dim_idx_2;
  for (i = 0; i < loop_ub; i++) {
    se->data[i] = 0.0;
  }

  for (i = 0; i < (int)dim_idx_0; i++) {
    for (loop_ub = 0; loop_ub < (int)dim_idx_1; loop_ub++) {
      for (k = 0; k < (int)dim_idx_2; k++) {
        dv5[0] = ((1.0 + (double)i) - 1.0) - r;
        dv5[1] = ((1.0 + (double)loop_ub) - 1.0) - r;
        dv5[2] = ((1.0 + (double)k) - 1.0) - r;
        if (norm(dv5) <= r) {
          se->data[(i + se->size[0] * loop_ub) + se->size[0] * se->size[1] * k] =
            1.0;
        }
      }
    }
  }
}
