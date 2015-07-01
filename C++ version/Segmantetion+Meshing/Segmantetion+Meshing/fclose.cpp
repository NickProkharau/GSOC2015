
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
#include "fclose.h"
#include "fileManager.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : double fileID
// Return Type  : void
//
void b_fclose(double fileID)
{
  c_fileManager(fileID);
}


