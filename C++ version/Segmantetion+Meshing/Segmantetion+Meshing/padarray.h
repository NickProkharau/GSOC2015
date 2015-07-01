
#ifndef __PADARRAY_H__
#define __PADARRAY_H__

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
extern void b_ConstantPad(const emxArray_uint8_T *a, signed char direction,
  emxArray_uint8_T *b);
extern void padarray(const emxArray_real_T *varargin_1, emxArray_real_T *b);

#endif


