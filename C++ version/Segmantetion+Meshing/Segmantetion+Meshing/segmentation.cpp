
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
#include "Segmentation1_emxutil.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void segmentation()
{
  emxArray_real_T *mri_anatomy;
  emxArray_real_T *segmentedmri;
  emxArray_real_T *b_segmentedmri;
  char expl_temp[8];
  char b_expl_temp[2];
  double c_expl_temp[16];
  double mri_dim[3];
  d_emxInit_real_T(&mri_anatomy, 7);
  c_emxInit_real_T(&segmentedmri, 3);
  c_emxInit_real_T(&b_segmentedmri, 3);

  //  reading MRI images
  b_read_mri_simple(mri_dim, mri_anatomy, c_expl_temp, b_expl_temp, expl_temp);
  //  segmenting MRI images
  b_volumesegment(mri_dim, mri_anatomy, segmentedmri);
  //  meshing MRI images
  b_meshing(segmentedmri, b_segmentedmri);
  emxFree_real_T(&b_segmentedmri);
  emxFree_real_T(&segmentedmri);
  emxFree_real_T(&mri_anatomy);
}


