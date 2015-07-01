
#ifndef __ROT_MATRIX_H__
#define __ROT_MATRIX_H__

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
extern void b_rot_matrix(const emxArray_real_T *Q, double M[16]);
extern void rot_matrix(double Q[3], double M[16]);

#endif
