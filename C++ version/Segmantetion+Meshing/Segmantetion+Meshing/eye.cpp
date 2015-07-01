
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
#include "eye.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : double I[12]
// Return Type  : void
//
void eye(double I[12])
{
  int k;
  memset(&I[0], 0, 12U * sizeof(double));
  for (k = 0; k < 3; k++) {
    I[k + (k << 2)] = 1.0;
  }
}
