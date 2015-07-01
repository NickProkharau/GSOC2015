
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
#include "prod.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : const double x[3]
// Return Type  : double
//
double prod(const double x[3])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 2; k++) {
    y *= x[k + 1];
  }

  return y;
}
