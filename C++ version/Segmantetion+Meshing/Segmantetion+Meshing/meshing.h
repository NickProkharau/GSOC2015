
#ifndef __MESHING_H__
#define __MESHING_H__

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
extern void b_meshing(const emxArray_real_T *vol, emxArray_real_T *b_vol);
extern void meshing(const boolean_T vol[8530021]);

#endif


