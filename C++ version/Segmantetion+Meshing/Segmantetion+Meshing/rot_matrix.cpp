
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
// Generate a rotation matrix from a quaternion xi+yj+zk+w,
//  where Q = [x y z], and w = 1-x^2-y^2-z^2.
// Arguments    : const emxArray_real_T *Q
//                double M[16]
// Return Type  : void
//
void b_rot_matrix(const emxArray_real_T *Q, double M[16])
{
  double y[3];
  int k;
  double x;
  double yw;
  double w;
  double z;
  double xx;
  double yy;
  double zz;
  double ww;
  double xy;
  double xz;
  double xw;
  double yz;
  static const signed char iv2[4] = { 0, 0, 0, 1 };

  //  Assume rigid body
  for (k = 0; k < 3; k++) {
    x = Q->data[k];
    y[k] = x * x;
  }

  yw = y[0];
  for (k = 0; k < 2; k++) {
    yw += y[k + 1];
  }

  w = sqrt(1.0 - yw);
  x = Q->data[0];
  yw = Q->data[1];
  z = Q->data[2];
  if (w < 1.0E-7) {
    w = 1.0 / sqrt((Q->data[0] * Q->data[0] + Q->data[1] * Q->data[1]) + Q->
                   data[2] * Q->data[2]);
    x = Q->data[0] * w;
    yw = Q->data[1] * w;
    z = Q->data[2] * w;
    w = 0.0;
  }

  xx = x * x;
  yy = yw * yw;
  zz = z * z;
  ww = w * w;
  xy = x * yw;
  xz = x * z;
  xw = x * w;
  yz = yw * z;
  yw *= w;
  x = z * w;
  M[0] = ((xx - yy) - zz) + ww;
  M[4] = 2.0 * (xy - x);
  M[8] = 2.0 * (xz + yw);
  M[12] = 0.0;
  M[1] = 2.0 * (xy + x);
  M[5] = ((-xx + yy) - zz) + ww;
  M[9] = 2.0 * (yz - xw);
  M[13] = 0.0;
  M[2] = 2.0 * (xz - yw);
  M[6] = 2.0 * (yz + xw);
  M[10] = ((-xx - yy) + zz) + ww;
  M[14] = 0.0;
  for (k = 0; k < 4; k++) {
    M[3 + (k << 2)] = iv2[k];
  }
}

//
// Generate a rotation matrix from a quaternion xi+yj+zk+w,
//  where Q = [x y z], and w = 1-x^2-y^2-z^2.
// Arguments    : double Q[3]
//                double M[16]
// Return Type  : void
//
void rot_matrix(double Q[3], double M[16])
{
  double y[3];
  int k;
  double yw;
  double w;
  double x;
  double z;
  double xx;
  double yy;
  double zz;
  double ww;
  double xy;
  double xz;
  double xw;
  double yz;
  static const signed char iv3[4] = { 0, 0, 0, 1 };

  //  Assume rigid body
  for (k = 0; k < 3; k++) {
    y[k] = Q[k] * Q[k];
  }

  yw = y[0];
  for (k = 0; k < 2; k++) {
    yw += y[k + 1];
  }

  w = sqrt(1.0 - yw);
  x = Q[0];
  yw = Q[1];
  z = Q[2];
  if (w < 1.0E-7) {
    w = 1.0 / sqrt((Q[0] * Q[0] + Q[1] * Q[1]) + Q[2] * Q[2]);
    x = Q[0] * w;
    yw = Q[1] * w;
    z = Q[2] * w;
    w = 0.0;
  }

  xx = x * x;
  yy = yw * yw;
  zz = z * z;
  ww = w * w;
  xy = x * yw;
  xz = x * z;
  xw = x * w;
  yz = yw * z;
  yw *= w;
  x = z * w;
  M[0] = ((xx - yy) - zz) + ww;
  M[4] = 2.0 * (xy - x);
  M[8] = 2.0 * (xz + yw);
  M[12] = 0.0;
  M[1] = 2.0 * (xy + x);
  M[5] = ((-xx + yy) - zz) + ww;
  M[9] = 2.0 * (yz - xw);
  M[13] = 0.0;
  M[2] = 2.0 * (xz - yw);
  M[6] = 2.0 * (yz + xw);
  M[10] = ((-xx - yy) + zz) + ww;
  M[14] = 0.0;
  for (k = 0; k < 4; k++) {
    M[3 + (k << 2)] = iv3[k];
  }
}

