
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
#include "diag.h"
#include "Segmentation1_emxutil.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : const emxArray_real_T *v
//                emxArray_real_T *d
// Return Type  : void
//
void diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int unnamed_idx_0;
  int unnamed_idx_1;
  int i21;
  unnamed_idx_0 = v->size[1];
  unnamed_idx_1 = v->size[1];
  i21 = d->size[0] * d->size[1];
  d->size[0] = unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)d, i21, (int)sizeof(double));
  i21 = d->size[0] * d->size[1];
  d->size[1] = unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)d, i21, (int)sizeof(double));
  unnamed_idx_0 *= unnamed_idx_1;
  for (i21 = 0; i21 < unnamed_idx_0; i21++) {
    d->data[i21] = 0.0;
  }

  for (unnamed_idx_0 = 0; unnamed_idx_0 + 1 <= v->size[1]; unnamed_idx_0++) {
    d->data[unnamed_idx_0 + d->size[0] * unnamed_idx_0] = v->data[unnamed_idx_0];
  }
}

