
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
#include "repmat.h"
#include "Segmentation1_emxutil.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : const int varargin_1[2]
//                emxArray_uint8_T *b
// Return Type  : void
//
void b_repmat(const int varargin_1[2], emxArray_uint8_T *b)
{
  int i38;
  int loop_ub;
  i38 = b->size[0] * b->size[1];
  b->size[0] = varargin_1[0];
  emxEnsureCapacity((emxArray__common *)b, i38, (int)sizeof(unsigned char));
  i38 = b->size[0] * b->size[1];
  b->size[1] = varargin_1[1];
  emxEnsureCapacity((emxArray__common *)b, i38, (int)sizeof(unsigned char));
  loop_ub = varargin_1[0] * varargin_1[1];
  for (i38 = 0; i38 < loop_ub; i38++) {
    b->data[i38] = 0;
  }
}

//
// Arguments    : const int varargin_1[2]
//                emxArray_real_T *b
// Return Type  : void
//
void repmat(const int varargin_1[2], emxArray_real_T *b)
{
  int i11;
  int loop_ub;
  i11 = b->size[0] * b->size[1];
  b->size[0] = varargin_1[0];
  emxEnsureCapacity((emxArray__common *)b, i11, (int)sizeof(double));
  i11 = b->size[0] * b->size[1];
  b->size[1] = varargin_1[1];
  emxEnsureCapacity((emxArray__common *)b, i11, (int)sizeof(double));
  loop_ub = varargin_1[0] * varargin_1[1];
  for (i11 = 0; i11 < loop_ub; i11++) {
    b->data[i11] = rtMinusInf;
  }
}

