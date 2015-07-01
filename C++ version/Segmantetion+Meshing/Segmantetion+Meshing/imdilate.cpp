
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
#include "imdilate.h"
#include "Segmentation1_emxutil.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : const emxArray_boolean_T *A
//                emxArray_boolean_T *B
// Return Type  : void
//
void imdilate(const emxArray_boolean_T *A, emxArray_boolean_T *B)
{
  int i34;
  int k;
  double asize[3];
  double nsize[3];
  static const boolean_T nhood[2197] = { false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, true, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, true,
    true, true, false, false, false, false, false, false, false, false, false,
    true, true, true, true, true, false, false, false, false, false, false,
    false, true, true, true, true, true, true, true, false, false, false, false,
    false, false, true, true, true, true, true, true, true, false, false, false,
    false, false, false, true, true, true, true, true, true, true, false, false,
    false, false, false, false, false, true, true, true, true, true, false,
    false, false, false, false, false, false, false, false, true, true, true,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, true, true, true,
    true, true, false, false, false, false, false, false, false, true, true,
    true, true, true, true, true, false, false, false, false, false, true, true,
    true, true, true, true, true, true, true, false, false, false, false, true,
    true, true, true, true, true, true, true, true, false, false, false, false,
    true, true, true, true, true, true, true, true, true, false, false, false,
    false, true, true, true, true, true, true, true, true, true, false, false,
    false, false, true, true, true, true, true, true, true, true, true, false,
    false, false, false, false, true, true, true, true, true, true, true, false,
    false, false, false, false, false, false, true, true, true, true, true,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, true, true, true, false, false, false, false,
    false, false, false, false, true, true, true, true, true, true, true, false,
    false, false, false, false, true, true, true, true, true, true, true, true,
    true, false, false, false, false, true, true, true, true, true, true, true,
    true, true, false, false, false, true, true, true, true, true, true, true,
    true, true, true, true, false, false, true, true, true, true, true, true,
    true, true, true, true, true, false, false, true, true, true, true, true,
    true, true, true, true, true, true, false, false, false, true, true, true,
    true, true, true, true, true, true, false, false, false, false, true, true,
    true, true, true, true, true, true, true, false, false, false, false, false,
    true, true, true, true, true, true, true, false, false, false, false, false,
    false, false, false, true, true, true, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, true, true, true,
    true, true, false, false, false, false, false, false, true, true, true, true,
    true, true, true, true, true, false, false, false, false, true, true, true,
    true, true, true, true, true, true, false, false, false, true, true, true,
    true, true, true, true, true, true, true, true, false, false, true, true,
    true, true, true, true, true, true, true, true, true, false, false, true,
    true, true, true, true, true, true, true, true, true, true, false, false,
    true, true, true, true, true, true, true, true, true, true, true, false,
    false, true, true, true, true, true, true, true, true, true, true, true,
    false, false, false, true, true, true, true, true, true, true, true, true,
    false, false, false, false, true, true, true, true, true, true, true, true,
    true, false, false, false, false, false, false, true, true, true, true, true,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    true, true, true, true, true, true, true, false, false, false, false, false,
    true, true, true, true, true, true, true, true, true, false, false, false,
    true, true, true, true, true, true, true, true, true, true, true, false,
    false, true, true, true, true, true, true, true, true, true, true, true,
    false, false, true, true, true, true, true, true, true, true, true, true,
    true, false, false, true, true, true, true, true, true, true, true, true,
    true, true, false, false, true, true, true, true, true, true, true, true,
    true, true, true, false, false, true, true, true, true, true, true, true,
    true, true, true, true, false, false, true, true, true, true, true, true,
    true, true, true, true, true, false, false, false, true, true, true, true,
    true, true, true, true, true, false, false, false, false, false, true, true,
    true, true, true, true, true, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, true, false, false, false, false, false,
    false, false, false, false, true, true, true, true, true, true, true, false,
    false, false, false, false, true, true, true, true, true, true, true, true,
    true, false, false, false, true, true, true, true, true, true, true, true,
    true, true, true, false, false, true, true, true, true, true, true, true,
    true, true, true, true, false, false, true, true, true, true, true, true,
    true, true, true, true, true, false, true, true, true, true, true, true,
    true, true, true, true, true, true, true, false, true, true, true, true,
    true, true, true, true, true, true, true, false, false, true, true, true,
    true, true, true, true, true, true, true, true, false, false, true, true,
    true, true, true, true, true, true, true, true, true, false, false, false,
    true, true, true, true, true, true, true, true, true, false, false, false,
    false, false, true, true, true, true, true, true, true, false, false, false,
    false, false, false, false, false, false, true, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, true, true, true, true,
    true, true, true, false, false, false, false, false, true, true, true, true,
    true, true, true, true, true, false, false, false, true, true, true, true,
    true, true, true, true, true, true, true, false, false, true, true, true,
    true, true, true, true, true, true, true, true, false, false, true, true,
    true, true, true, true, true, true, true, true, true, false, false, true,
    true, true, true, true, true, true, true, true, true, true, false, false,
    true, true, true, true, true, true, true, true, true, true, true, false,
    false, true, true, true, true, true, true, true, true, true, true, true,
    false, false, true, true, true, true, true, true, true, true, true, true,
    true, false, false, false, true, true, true, true, true, true, true, true,
    true, false, false, false, false, false, true, true, true, true, true, true,
    true, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, true, true, true, true, true, false, false, false, false, false,
    false, true, true, true, true, true, true, true, true, true, false, false,
    false, false, true, true, true, true, true, true, true, true, true, false,
    false, false, true, true, true, true, true, true, true, true, true, true,
    true, false, false, true, true, true, true, true, true, true, true, true,
    true, true, false, false, true, true, true, true, true, true, true, true,
    true, true, true, false, false, true, true, true, true, true, true, true,
    true, true, true, true, false, false, true, true, true, true, true, true,
    true, true, true, true, true, false, false, false, true, true, true, true,
    true, true, true, true, true, false, false, false, false, true, true, true,
    true, true, true, true, true, true, false, false, false, false, false, false,
    true, true, true, true, true, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, true, true, true, false,
    false, false, false, false, false, false, false, true, true, true, true,
    true, true, true, false, false, false, false, false, true, true, true, true,
    true, true, true, true, true, false, false, false, false, true, true, true,
    true, true, true, true, true, true, false, false, false, true, true, true,
    true, true, true, true, true, true, true, true, false, false, true, true,
    true, true, true, true, true, true, true, true, true, false, false, true,
    true, true, true, true, true, true, true, true, true, true, false, false,
    false, true, true, true, true, true, true, true, true, true, false, false,
    false, false, true, true, true, true, true, true, true, true, true, false,
    false, false, false, false, true, true, true, true, true, true, true, false,
    false, false, false, false, false, false, false, true, true, true, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, true, true, true, true, true, false, false, false,
    false, false, false, false, true, true, true, true, true, true, true, false,
    false, false, false, false, true, true, true, true, true, true, true, true,
    true, false, false, false, false, true, true, true, true, true, true, true,
    true, true, false, false, false, false, true, true, true, true, true, true,
    true, true, true, false, false, false, false, true, true, true, true, true,
    true, true, true, true, false, false, false, false, true, true, true, true,
    true, true, true, true, true, false, false, false, false, false, true, true,
    true, true, true, true, true, false, false, false, false, false, false,
    false, true, true, true, true, true, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, true, true, true, false, false, false, false, false,
    false, false, false, false, true, true, true, true, true, false, false,
    false, false, false, false, false, true, true, true, true, true, true, true,
    false, false, false, false, false, false, true, true, true, true, true, true,
    true, false, false, false, false, false, false, true, true, true, true, true,
    true, true, false, false, false, false, false, false, false, true, true,
    true, true, true, false, false, false, false, false, false, false, false,
    false, true, true, true, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    true, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false };

  for (i34 = 0; i34 < 3; i34++) {
    k = B->size[0] * B->size[1] * B->size[2];
    B->size[i34] = A->size[i34];
    emxEnsureCapacity((emxArray__common *)B, k, (int)sizeof(boolean_T));
  }

  k = 3;
  while ((k > 2) && (A->size[2] == 1)) {
    k = 2;
  }

  for (i34 = 0; i34 < 3; i34++) {
    asize[i34] = A->size[i34];
  }

  for (i34 = 0; i34 < 3; i34++) {
    nsize[i34] = 13.0;
  }

  //dilate_binary_tbb(&A->data[0], asize, (double)k, nhood, nsize, 3.0, &B->data[0]);
}

