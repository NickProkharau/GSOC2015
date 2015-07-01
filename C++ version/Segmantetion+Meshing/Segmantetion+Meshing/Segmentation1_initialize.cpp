
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
#include "Segmentation1_initialize.h"
#include "fileManager.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void Segmentation1_initialize()
{
  rt_InitInfAndNaN(8U);
  filedata_init();
}

