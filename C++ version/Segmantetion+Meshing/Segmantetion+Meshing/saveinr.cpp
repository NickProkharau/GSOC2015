
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
#include "fopen.h"
#include <stdio.h>


// Function Definitions

//
// saveinr(vol,fname)
//
//  save a surface mesh to INR Format
//
//  input:
//       vol: input, a binary volume
//       fname: output file name
// Arguments    : const unsigned char vol[8530021]
//                const emxArray_char_T *fname
// Return Type  : void
//
void saveinr(const unsigned char [8530021], const emxArray_char_T *fname)
{
  double fid;
  fid = h_fopen(fname);

  // fwrite(fid,header,'char');
  // fwrite(fid,vol,dtype);
  b_fclose(fid);
}

