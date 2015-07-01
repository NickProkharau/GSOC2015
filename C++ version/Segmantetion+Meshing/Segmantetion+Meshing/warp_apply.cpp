
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
// FT_WARP_APPLY performs a 3D linear or nonlinear transformation on the input
//  coordinates
// Arguments    : const double M[16]
//                const double input[24]
//                const emxArray_char_T *method
//                double tol
//                double warped[24]
// Return Type  : void
//
void warp_apply(const double M[16], const double input[24], const
                emxArray_char_T *method, double tol, double warped[24])
{
  boolean_T b_bool;
  int k;
  int32_T exitg10;
  static const char cv48[9] = { 'n', 'o', 'n', 'l', 'i', 'n', 'e', 'a', 'r' };

  boolean_T c_bool;
  int32_T exitg9;
  static const char cv49[7] = { 'n', 'o', 'n', 'l', 'i', 'n', '1' };

  boolean_T d_bool;
  int32_T exitg8;
  static const char cv50[7] = { 'n', 'o', 'n', 'l', 'i', 'n', '0' };

  boolean_T e_bool;
  int32_T exitg7;
  static const char cv51[7] = { 'n', 'o', 'n', 'l', 'i', 'n', '2' };

  boolean_T f_bool;
  int32_T exitg6;
  static const char cv52[7] = { 'n', 'o', 'n', 'l', 'i', 'n', '3' };

  boolean_T g_bool;
  int32_T exitg5;
  static const char cv53[7] = { 'n', 'o', 'n', 'l', 'i', 'n', '4' };

  boolean_T h_bool;
  int32_T exitg4;
  static const char cv54[7] = { 'n', 'o', 'n', 'l', 'i', 'n', '5' };

  boolean_T x[7];
  boolean_T exitg3;
  int32_T exitg2;
  static const char cv55[10] = { 'h', 'o', 'm', 'o', 'g', 'e', 'n', 'o', 'u',
    's' };

  boolean_T guard1 = false;
  int32_T exitg1;
  static const char cv56[11] = { 'h', 'o', 'm', 'o', 'g', 'e', 'n', 'e', 'o',
    'u', 's' };

  double b_input[32];
  int i40;
  int i41;
  double b_warped;
  double c_warped;

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  //  nonlinear warping
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  b_bool = false;
  if (method->size[1] != 9) {
  } else {
    k = 0;
    do {
      exitg10 = 0;
      if (k <= 8) {
        if (method->data[k] != cv48[k]) {
          exitg10 = 1;
        } else {
          k++;
        }
      } else {
        b_bool = true;
        exitg10 = 1;
      }
    } while (exitg10 == 0);
  }

  c_bool = false;
  if (method->size[1] != 7) {
  } else {
    k = 0;
    do {
      exitg9 = 0;
      if (k <= 6) {
        if (method->data[k] != cv49[k]) {
          exitg9 = 1;
        } else {
          k++;
        }
      } else {
        c_bool = true;
        exitg9 = 1;
      }
    } while (exitg9 == 0);
  }

  d_bool = false;
  if (method->size[1] != 7) {
  } else {
    k = 0;
    do {
      exitg8 = 0;
      if (k <= 6) {
        if (method->data[k] != cv50[k]) {
          exitg8 = 1;
        } else {
          k++;
        }
      } else {
        d_bool = true;
        exitg8 = 1;
      }
    } while (exitg8 == 0);
  }

  e_bool = false;
  if (method->size[1] != 7) {
  } else {
    k = 0;
    do {
      exitg7 = 0;
      if (k <= 6) {
        if (method->data[k] != cv51[k]) {
          exitg7 = 1;
        } else {
          k++;
        }
      } else {
        e_bool = true;
        exitg7 = 1;
      }
    } while (exitg7 == 0);
  }

  f_bool = false;
  if (method->size[1] != 7) {
  } else {
    k = 0;
    do {
      exitg6 = 0;
      if (k <= 6) {
        if (method->data[k] != cv52[k]) {
          exitg6 = 1;
        } else {
          k++;
        }
      } else {
        f_bool = true;
        exitg6 = 1;
      }
    } while (exitg6 == 0);
  }

  g_bool = false;
  if (method->size[1] != 7) {
  } else {
    k = 0;
    do {
      exitg5 = 0;
      if (k <= 6) {
        if (method->data[k] != cv53[k]) {
          exitg5 = 1;
        } else {
          k++;
        }
      } else {
        g_bool = true;
        exitg5 = 1;
      }
    } while (exitg5 == 0);
  }

  h_bool = false;
  if (method->size[1] != 7) {
  } else {
    k = 0;
    do {
      exitg4 = 0;
      if (k <= 6) {
        if (method->data[k] != cv54[k]) {
          exitg4 = 1;
        } else {
          k++;
        }
      } else {
        h_bool = true;
        exitg4 = 1;
      }
    } while (exitg4 == 0);
  }

  x[0] = b_bool;
  x[1] = c_bool;
  x[2] = d_bool;
  x[3] = e_bool;
  x[4] = f_bool;
  x[5] = g_bool;
  x[6] = h_bool;
  b_bool = false;
  k = 0;
  exitg3 = false;
  while ((!exitg3) && (k < 7)) {
    if (!(x[k] == 0)) {
      b_bool = true;
      exitg3 = true;
    } else {
      k++;
    }
  }

  if (b_bool) {
  } else {
    b_bool = false;
    if (method->size[1] != 10) {
    } else {
      k = 0;
      do {
        exitg2 = 0;
        if (k <= 9) {
          if (method->data[k] != cv55[k]) {
            exitg2 = 1;
          } else {
            k++;
          }
        } else {
          b_bool = true;
          exitg2 = 1;
        }
      } while (exitg2 == 0);
    }

    guard1 = false;
    if (b_bool) {
      guard1 = true;
    } else {
      b_bool = false;
      if (method->size[1] != 11) {
      } else {
        k = 0;
        do {
          exitg1 = 0;
          if (k <= 10) {
            if (method->data[k] != cv56[k]) {
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
        guard1 = true;
      }
    }

    if (guard1) {
      // warped = M * [input'; ones(1, size(input, 1))];
      // warped = warped(1:3,:)';
      //  below achieves the same as lines 154-155
      for (k = 0; k < 3; k++) {
        memcpy(&b_input[k << 3], &input[k << 3], sizeof(double) << 3);
      }

      for (k = 0; k < 8; k++) {
        b_input[24 + k] = 1.0;
        for (i40 = 0; i40 < 3; i40++) {
          warped[k + (i40 << 3)] = 0.0;
          for (i41 = 0; i41 < 4; i41++) {
            warped[k + (i40 << 3)] += b_input[k + (i41 << 3)] * M[i40 + (i41 <<
              2)];
          }
        }
      }

      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      //  using external function that returns a homogeneous transformation matrix 
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      //  elseif exist(method, 'file') && ~isa(M, 'struct')
      //    % get the homogenous transformation matrix
      //    H = feval(method, M);
      //    warped = ft_warp_apply(H, input, 'homogeneous');
    }
  }

  if (tol > 0.0) {
    for (k = 0; k < 24; k++) {
      b_warped = warped[k] / tol;
      if (b_warped < 0.0) {
        c_warped = ceil(b_warped);
      } else {
        c_warped = floor(b_warped);
      }

      warped[k] = c_warped * tol;
    }
  }
}

