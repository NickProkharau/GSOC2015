
#ifndef __SEGMENTATION1_TYPES_H__
#define __SEGMENTATION1_TYPES_H__

/* Include files */
#include "stdafx.h"
#include "rtwtypes.h"

/* Type Definitions */
#include <stdio.h>
#ifndef struct_emxArray__common
#define struct_emxArray__common
struct emxArray__common
{
    void *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray__common*/
#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T
struct emxArray_boolean_T
{
    boolean_T *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_boolean_T*/
#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T
struct emxArray_char_T
{
    char *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_char_T*/
#ifndef struct_emxArray_int16_T
#define struct_emxArray_int16_T
struct emxArray_int16_T
{
    short *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_int16_T*/
#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T
struct emxArray_int32_T
{
    int *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_int32_T*/
#ifndef struct_emxArray_real32_T
#define struct_emxArray_real32_T
struct emxArray_real32_T
{
    float *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_real32_T*/
#ifndef struct_emxArray_uint8_T
#define struct_emxArray_uint8_T
struct emxArray_uint8_T
{
    unsigned char *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_uint8_T*/
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T
{
    double *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_real_T*/
typedef struct
{
    double dim[3];
    emxArray_real_T *anatomy;
    double transform[16];
    char unit[2];
    char coordsys[8];
} struct0_T;
typedef struct
{
    double dim[3];
    double anatomy[8530021];
    double transform[16];
    char unit[2];
    char coordsys[8];
} struct1_T;

#endif
