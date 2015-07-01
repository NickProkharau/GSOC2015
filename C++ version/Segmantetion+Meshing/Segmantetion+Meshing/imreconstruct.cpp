
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
#include "imreconstruct.h"
#include <stdio.h>

// Function Definitions

//
// Arguments    : emxArray_real_T *varargin_1
//                const emxArray_real_T *varargin_2
// Return Type  : void
//
void imreconstruct(emxArray_real_T *varargin_1, const emxArray_real_T
                   *varargin_2)
{
  double imSize[2];
  int i42;
  boolean_T conn[9];
  double connSize[2];
  for (i42 = 0; i42 < 2; i42++) {
    imSize[i42] = varargin_1->size[i42];
  }

  for (i42 = 0; i42 < 9; i42++) {
    conn[i42] = true;
  }

  for (i42 = 0; i42 < 2; i42++) {
    connSize[i42] = 3.0;
  }

  //imreconstruct_real64(&varargin_1->data[0], &varargin_2->data[0], 2.0, imSize,
                       //conn, 2.0, connSize);
}


