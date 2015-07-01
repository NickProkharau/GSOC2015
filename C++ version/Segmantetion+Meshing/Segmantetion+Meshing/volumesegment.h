
#ifndef __VOLUMESEGMENT_H__
#define __VOLUMESEGMENT_H__

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
extern void b_volumesegment(const double mri_dim[3], const emxArray_real_T
  *mri_anatomy, emxArray_real_T *segmented);
extern void volumesegment(const struct1_T *mri, double segmented[8530021]);

#endif

