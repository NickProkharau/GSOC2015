
// Include files
#include "stdafx.h"
#include "rt_nonfinite.h"
#include "cgalv2m.h"
#include "corner_points.h"
#include "fillholes.h"
#include "fillimg.h"
#include "meshing.h"
#include "read_mri_simple.h"
#include "readmedit.h"
#include "rot_matrix.h"
#include "saveinr.h"
#include "segmentation.h"
#include "sphere.h"
#include "volumesegment.h"
#include "warp_apply.h"
#include "fread.h"
#include "Segmentation1_emxutil.h"
#include "fileManager.h"
#include <stdio.h>


// Function Definitions

//
// Arguments    : double fileID
//                emxArray_real_T *A
// Return Type  : void
//
void b_fread(double fileID, emxArray_real_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_real_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_real_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
  } else {
    numRead = A->size[0];
    A->size[0] = 100564;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(double));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 100564) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 8, 100564 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 100565; numRead++) {
      A->data[numRead] = 0.0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 100564) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(double));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_real_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_real_T *A
// Return Type  : void
//
void c_fread(double fileID, emxArray_real_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_real_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_real_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
  } else {
    numRead = A->size[0];
    A->size[0] = 565750;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(double));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 565750) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 8, 565750 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 565751; numRead++) {
      A->data[numRead] = 0.0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 565750) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(double));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_real_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_real_T *A
// Return Type  : void
//
void d_fread(double fileID, emxArray_real_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_real_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_real_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
  } else {
    numRead = A->size[0];
    A->size[0] = 362368;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(double));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 362368) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 8, 362368 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 362369; numRead++) {
      A->data[numRead] = 0.0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 362368) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(double));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(double));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_real_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_int32_T *A
// Return Type  : void
//
void e_fread(double fileID, emxArray_int32_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_int32_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  b_emxInit_int32_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(int));
  } else {
    numRead = A->size[0];
    A->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(int));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 1) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 4, 1 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    numRead = bytesOut + 1;
    while (numRead <= 1) {
      A->data[0] = 0;
      numRead = 2;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(int));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 1) {
      numRead = b_A->size[0];
      b_A->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(int));
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(int));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_int32_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_uint8_T *A
// Return Type  : void
//
void f_fread(double fileID, emxArray_uint8_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_uint8_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_uint8_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
  } else {
    numRead = A->size[0];
    A->size[0] = 10;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 10) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 1, 10 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 11; numRead++) {
      A->data[numRead] = 0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 10) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
        char));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_uint8_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_uint8_T *A
// Return Type  : void
//
void g_fread(double fileID, emxArray_uint8_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_uint8_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_uint8_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
  } else {
    numRead = A->size[0];
    A->size[0] = 18;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 18) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 1, 18 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 19; numRead++) {
      A->data[numRead] = 0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 18) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
        char));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_uint8_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_int16_T *A
// Return Type  : void
//
void h_fread(double fileID, emxArray_int16_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_int16_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  b_emxInit_int16_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(short));
  } else {
    numRead = A->size[0];
    A->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(short));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 1) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 2, 1 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    numRead = bytesOut + 1;
    while (numRead <= 1) {
      A->data[0] = 0;
      numRead = 2;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(short));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 1) {
      numRead = b_A->size[0];
      b_A->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(short));
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(short));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_int16_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_uint8_T *A
// Return Type  : void
//
void i_fread(double fileID, emxArray_uint8_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_uint8_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_uint8_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
  } else {
    numRead = A->size[0];
    A->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 1) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 1, 1 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    numRead = bytesOut + 1;
    while (numRead <= 1) {
      A->data[0] = 0;
      numRead = 2;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 1) {
      numRead = b_A->size[0];
      b_A->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
        char));
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_uint8_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_int16_T *A
// Return Type  : void
//
void j_fread(double fileID, emxArray_int16_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_int16_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  b_emxInit_int16_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(short));
  } else {
    numRead = A->size[0];
    A->size[0] = 8;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(short));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 8) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 2, 8 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 9; numRead++) {
      A->data[numRead] = 0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(short));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 8) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(short));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(short));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_int16_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_real32_T *A
// Return Type  : void
//
void k_fread(double fileID, emxArray_real32_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_real32_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  b_emxInit_real32_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
  } else {
    numRead = A->size[0];
    A->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(float));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 1) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 4, 1 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    numRead = bytesOut + 1;
    while (numRead <= 1) {
      A->data[0] = 0.0F;
      numRead = 2;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 1) {
      numRead = b_A->size[0];
      b_A->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(float));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_real32_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_real32_T *A
// Return Type  : void
//
void l_fread(double fileID, emxArray_real32_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_real32_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  b_emxInit_real32_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
  } else {
    numRead = A->size[0];
    A->size[0] = 8;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(float));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 8) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 4, 8 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 9; numRead++) {
      A->data[numRead] = 0.0F;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 8) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(float));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_real32_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_uint8_T *A
// Return Type  : void
//
void m_fread(double fileID, emxArray_uint8_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_uint8_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_uint8_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
  } else {
    numRead = A->size[0];
    A->size[0] = 80;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 80) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 1, 80 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 81; numRead++) {
      A->data[numRead] = 0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 80) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
        char));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_uint8_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_uint8_T *A
// Return Type  : void
//
void n_fread(double fileID, emxArray_uint8_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_uint8_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_uint8_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
  } else {
    numRead = A->size[0];
    A->size[0] = 24;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 24) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 1, 24 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 25; numRead++) {
      A->data[numRead] = 0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 24) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
        char));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_uint8_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_real32_T *A
// Return Type  : void
//
void o_fread(double fileID, emxArray_real32_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_real32_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  b_emxInit_real32_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
  } else {
    numRead = A->size[0];
    A->size[0] = 4;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(float));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 4) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 4, 4 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 5; numRead++) {
      A->data[numRead] = 0.0F;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 4) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(float));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(float));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_real32_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_uint8_T *A
// Return Type  : void
//
void p_fread(double fileID, emxArray_uint8_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_uint8_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_uint8_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
  } else {
    numRead = A->size[0];
    A->size[0] = 16;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 16) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 1, 16 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 17; numRead++) {
      A->data[numRead] = 0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 16) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
        char));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_uint8_T(&b_A);
}

//
// Arguments    : double fileID
//                emxArray_uint8_T *A
// Return Type  : void
//
void q_fread(double fileID, emxArray_uint8_T *A)
{
  FILE * filestar;
  boolean_T p;
  emxArray_uint8_T *b_A;
  int numRead;
  int bytesOut;
  size_t numReadSizeT;
  int loop_ub;
  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_uint8_T(&b_A, 1);
  if (filestar == (FILE *)NULL) {
    numRead = b_A->size[0];
    b_A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
  } else {
    numRead = A->size[0];
    A->size[0] = 4;
    emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
    bytesOut = 0;
    numRead = 1;
    while ((bytesOut < 4) && (numRead > 0)) {
      numReadSizeT = fread(&A->data[0], 1, 4 - bytesOut, filestar);
      numRead = (int)numReadSizeT;
      bytesOut += (int)numReadSizeT;
    }

    for (numRead = bytesOut; numRead + 1 < 5; numRead++) {
      A->data[numRead] = 0;
    }

    numRead = b_A->size[0];
    b_A->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
      char));
    loop_ub = A->size[0];
    for (numRead = 0; numRead < loop_ub; numRead++) {
      b_A->data[numRead] = A->data[numRead];
    }

    if (bytesOut < 4) {
      if (1 > bytesOut) {
        loop_ub = 0;
      } else {
        loop_ub = bytesOut;
      }

      numRead = b_A->size[0];
      b_A->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)b_A, numRead, (int)sizeof(unsigned
        char));
      for (numRead = 0; numRead < loop_ub; numRead++) {
        b_A->data[numRead] = A->data[numRead];
      }
    }
  }

  numRead = A->size[0];
  A->size[0] = b_A->size[0];
  emxEnsureCapacity((emxArray__common *)A, numRead, (int)sizeof(unsigned char));
  loop_ub = b_A->size[0];
  for (numRead = 0; numRead < loop_ub; numRead++) {
    A->data[numRead] = b_A->data[numRead];
  }

  emxFree_uint8_T(&b_A);
}

//
// Arguments    : double fileID
//                double sizeA
//                emxArray_real_T *A
// Return Type  : void
//
void r_fread(double fileID, double sizeA, emxArray_real_T *A)
{
  int dims_idx_0;
  boolean_T doEOF;
  FILE * filestar;
  boolean_T p;
  emxArray_real_T *b_A;
  int i22;
  int c;
  int numRead;
  size_t numReadSizeT;
  double tbuf[1024];
  if (rtIsInf(sizeA)) {
    dims_idx_0 = 1024;
    doEOF = true;
  } else {
    dims_idx_0 = (int)sizeA;
    doEOF = false;
  }

  filestar = f_fileManager(fileID);
  if ((fileID != 0.0) && (fileID != 1.0) && (fileID != 2.0)) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    filestar = NULL;
  }

  emxInit_real_T(&b_A, 1);
  if (!doEOF) {
    if (filestar == (FILE *)NULL) {
      i22 = A->size[0];
      A->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)A, i22, (int)sizeof(double));
    } else {
      i22 = b_A->size[0];
      b_A->size[0] = (int)sizeA;
      emxEnsureCapacity((emxArray__common *)b_A, i22, (int)sizeof(double));
      c = 0;
      numRead = 1;
      while ((c < dims_idx_0) && (numRead > 0)) {
        numReadSizeT = fread(&b_A->data[0], 8, dims_idx_0 - c, filestar);
        numRead = (int)numReadSizeT;
        c += (int)numReadSizeT;
      }

      i22 = b_A->size[0];
      for (dims_idx_0 = c; dims_idx_0 + 1 <= i22; dims_idx_0++) {
        b_A->data[dims_idx_0] = 0.0;
      }

      i22 = A->size[0];
      A->size[0] = b_A->size[0];
      emxEnsureCapacity((emxArray__common *)A, i22, (int)sizeof(double));
      dims_idx_0 = b_A->size[0];
      for (i22 = 0; i22 < dims_idx_0; i22++) {
        A->data[i22] = b_A->data[i22];
      }

      if (c < sizeA) {
        if (1 > c) {
          dims_idx_0 = 0;
        } else {
          dims_idx_0 = c;
        }

        i22 = A->size[0];
        A->size[0] = dims_idx_0;
        emxEnsureCapacity((emxArray__common *)A, i22, (int)sizeof(double));
        for (i22 = 0; i22 < dims_idx_0; i22++) {
          A->data[i22] = b_A->data[i22];
        }
      }
    }
  } else {
    i22 = A->size[0];
    A->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)A, i22, (int)sizeof(double));
    if (filestar == (FILE *)NULL) {
    } else {
      c = 1;
      while (c > 0) {
        c = 0;
        numRead = 1;
        while ((c < 1024) && (numRead > 0)) {
          numReadSizeT = fread(tbuf, 8, 1024 - c, filestar);
          numRead = (int)numReadSizeT;
          c += (int)numReadSizeT;
        }

        if (1 > c) {
          dims_idx_0 = 0;
        } else {
          dims_idx_0 = c;
        }

        numRead = A->size[0];
        i22 = A->size[0];
        A->size[0] = numRead + dims_idx_0;
        emxEnsureCapacity((emxArray__common *)A, i22, (int)sizeof(double));
        for (i22 = 0; i22 < dims_idx_0; i22++) {
          A->data[numRead + i22] = tbuf[i22];
        }
      }
    }
  }

  emxFree_real_T(&b_A);
}
