
#ifndef __FILLHOLES_H__
#define __FILLHOLES_H__

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
extern void b_fillholes(const emxArray_real_T *input, emxArray_boolean_T *output);
extern void c_fillholes(const emxArray_boolean_T *input, emxArray_boolean_T
  *output);
extern void fillholes(const double input[8530021], boolean_T output[8530021]);

#endif

