
#ifndef __READ_MRI_SIMPLE_H__
#define __READ_MRI_SIMPLE_H__

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
extern void b_read_mri_simple(double mri_dim[3], emxArray_real_T *mri_anatomy,
  double mri_transform[16], char mri_unit[2], char mri_coordsys[8]);
extern void read_mri_simple(const emxArray_char_T *fname, struct0_T *mri);

#endif
