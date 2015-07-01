
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
#include "fopen.h"
#include "fileManager.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : void
// Return Type  : double
//
double b_fopen()
{
  return b_fileManager();
}

//
// Arguments    : void
// Return Type  : double
//
double c_fopen()
{
  return e_fileManager();
}

//
// Arguments    : void
// Return Type  : double
//
double d_fopen()
{
  return g_fileManager();
}

//
// Arguments    : void
// Return Type  : double
//
double e_fopen()
{
  return h_fileManager();
}

//
// Arguments    : const emxArray_char_T *filename
// Return Type  : double
//
double f_fopen(const emxArray_char_T *filename)
{
  double fileID;
  boolean_T b_bool;
  int k;
  int32_T exitg1;
  static const char cv19[3] = { 'a', 'l', 'l' };

  b_bool = false;
  if (filename->size[1] != 3) {
  } else {
    k = 0;
    do {
      exitg1 = 0;
      if (k <= 2) {
        if (filename->data[k] != cv19[k]) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    fileID = 0.0;
  } else {
    fileID = i_fileManager(filename);
  }

  return fileID;
}

//
// Arguments    : const emxArray_char_T *filename
// Return Type  : double
//
double g_fopen(const emxArray_char_T *filename)
{
  return i_fileManager(filename);
}

//
// Arguments    : const emxArray_char_T *filename
// Return Type  : double
//
double h_fopen(const emxArray_char_T *filename)
{
  return j_fileManager(filename);
}

//
// Arguments    : void
// Return Type  : double
//
double i_fopen()
{
  return k_fileManager();
}

//
// Arguments    : void
// Return Type  : double
//
double j_fopen()
{
  return l_fileManager();
}

//
// Arguments    : void
// Return Type  : double
//
double k_fopen()
{
  return m_fileManager();
}

//
// Arguments    : void
// Return Type  : double
//
double l_fopen()
{
  return n_fileManager();
}

//
// Arguments    : void
// Return Type  : double
//
double m_fopen()
{
  return o_fileManager();
}

