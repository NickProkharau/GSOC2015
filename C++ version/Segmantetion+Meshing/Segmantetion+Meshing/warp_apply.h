
#ifndef __WARP_APPLY_H__
#define __WARP_APPLY_H__

// Include files
#include "stdafx.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "Segmentation1_types.h"

// Function Declarations
extern void warp_apply(const double M[16], const double input[24], const
  emxArray_char_T *method, double tol, double warped[24]);

#endif

