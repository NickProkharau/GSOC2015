
#ifndef __FREAD_H__
#define __FREAD_H__

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
extern void b_fread(double fileID, emxArray_real_T *A);
extern void c_fread(double fileID, emxArray_real_T *A);
extern void d_fread(double fileID, emxArray_real_T *A);
extern void e_fread(double fileID, emxArray_int32_T *A);
extern void f_fread(double fileID, emxArray_uint8_T *A);
extern void g_fread(double fileID, emxArray_uint8_T *A);
extern void h_fread(double fileID, emxArray_int16_T *A);
extern void i_fread(double fileID, emxArray_uint8_T *A);
extern void j_fread(double fileID, emxArray_int16_T *A);
extern void k_fread(double fileID, emxArray_real32_T *A);
extern void l_fread(double fileID, emxArray_real32_T *A);
extern void m_fread(double fileID, emxArray_uint8_T *A);
extern void n_fread(double fileID, emxArray_uint8_T *A);
extern void o_fread(double fileID, emxArray_real32_T *A);
extern void p_fread(double fileID, emxArray_uint8_T *A);
extern void q_fread(double fileID, emxArray_uint8_T *A);
extern void r_fread(double fileID, double sizeA, emxArray_real_T *A);

#endif

