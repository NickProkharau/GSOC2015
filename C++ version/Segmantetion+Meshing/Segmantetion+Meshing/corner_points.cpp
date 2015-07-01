
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
#include <stdio.h>

// Function Definitions

//
// corner_points returns the eight corner points of an anatomical volume
//  in voxel and in head coordinates
// Arguments    : const double dim[3]
//                const double transform[16]
//                double head[24]
// Return Type  : void
//
void corner_points(const double dim[3], const double transform[16], double head
                   [24])
{
  double dv0[32];
  int i4;
  int i5;
  int i6;

  //  determine the corner points of the volume in voxel space
  //  determine the corner points of the volume in plotting space
  //  FT_WARP_APPLY performs a 3D linear or nonlinear transformation on the input 
  //  coordinates
  //  no specific transformation mode has been selected
  //  it looks like a homogenous transformation matrix
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  //  nonlinear warping
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  //  below achieves the same as lines 154-155
  for (i4 = 0; i4 < 3; i4++) {
    dv0[i4 << 3] = 1.0;
  }

  dv0[1] = dim[0];
  dv0[9] = 1.0;
  dv0[17] = 1.0;
  dv0[2] = dim[0];
  dv0[10] = dim[1];
  dv0[18] = 1.0;
  dv0[3] = 1.0;
  dv0[11] = dim[1];
  dv0[19] = 1.0;
  dv0[4] = 1.0;
  dv0[12] = 1.0;
  dv0[20] = dim[2];
  dv0[5] = dim[0];
  dv0[13] = 1.0;
  dv0[21] = dim[2];
  dv0[6] = dim[0];
  dv0[14] = dim[1];
  dv0[22] = dim[2];
  dv0[7] = 1.0;
  dv0[15] = dim[1];
  dv0[23] = dim[2];
  for (i4 = 0; i4 < 8; i4++) {
    dv0[24 + i4] = 1.0;
    for (i5 = 0; i5 < 3; i5++) {
      head[i4 + (i5 << 3)] = 0.0;
      for (i6 = 0; i6 < 4; i6++) {
        head[i4 + (i5 << 3)] += dv0[i4 + (i6 << 3)] * transform[i5 + (i6 << 2)];
      }
    }
  }
}


