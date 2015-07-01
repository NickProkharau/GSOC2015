
#ifndef __FILLIMG_H__
#define __FILLIMG_H__

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
extern void b_fillimg(const double I[44197], emxArray_real_T *I2);
extern void c_fillimg(const double I[37249], emxArray_real_T *I2);
extern void d_fillimg(const double I[44197], emxArray_real_T *I2);
extern void e_fillimg(const emxArray_real_T *I, emxArray_real_T *I2);
extern void f_fillimg(const emxArray_boolean_T *I, emxArray_boolean_T *I2);
extern void fillimg(const emxArray_real_T *I, emxArray_real_T *I2);

#endif

