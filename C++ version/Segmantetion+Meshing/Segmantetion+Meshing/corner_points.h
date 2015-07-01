
#ifndef __CORNER_POINTS_H__
#define __CORNER_POINTS_H__

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
extern void corner_points(const double dim[3], const double transform[16],
  double head[24]);

#endif

