
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
#include "fileManager.h"
#include "Segmentation1_emxutil.h"
#include "Segmentation1_rtwutil.h"
#include <stdio.h>


// Variable Definitions
static FILE * eml_openfiles[20];
static boolean_T eml_autoflush[20];

// Function Declarations
static FILE * d_fileManager(signed char varargin_1);
static signed char filedata();

// Function Definitions

//
// Arguments    : signed char varargin_1
// Return Type  : FILE *
//
static FILE * d_fileManager(signed char varargin_1)
{
  FILE * f;
  signed char fileid;
  fileid = varargin_1;
  if ((varargin_1 > 22) || (varargin_1 < 0)) {
    fileid = -1;
  }

  if (fileid >= 3) {
    f = eml_openfiles[fileid - 3];
  } else if (fileid == 0) {
    f = stdin;
  } else if (fileid == 1) {
    f = stdout;
  } else if (fileid == 2) {
    f = stderr;
  } else {
    f = NULL;
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
static signed char filedata()
{
  signed char f;
  signed char k;
  boolean_T exitg1;
  f = 0;
  k = 1;
  exitg1 = false;
  while ((!exitg1) && (k < 21)) {
    if (eml_openfiles[k - 1] == (FILE *)NULL) {
      f = k;
      exitg1 = true;
    } else {
      k++;
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char b_fileManager()
{
  signed char f;
  signed char j;
  char cv0[17];
  int i0;
  static const char cv1[17] = { 'p', 'r', 'e', '_', 'c', 'g', 'a', 'l', 'm', 'e',
    's', 'h', '.', 'i', 'n', 'r', '\x00' };

  char cv2[3];
  static const char cv3[3] = { 'w', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i0 = 0; i0 < 17; i0++) {
      cv0[i0] = cv1[i0];
    }

    for (i0 = 0; i0 < 3; i0++) {
      cv2[i0] = cv3[i0];
    }

    filestar = fopen(cv0, cv2);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : double varargin_1
// Return Type  : int
//
int c_fileManager(double varargin_1)
{
  int f;
  signed char fileid;
  FILE * filestar;
  int cst;
  f = -1;
  fileid = (signed char)rt_roundd_snf(varargin_1);
  if ((fileid > 22) || (fileid < 0) || (varargin_1 != fileid)) {
    fileid = -1;
  }

  filestar = d_fileManager(fileid);
  if ((filestar == (FILE *)NULL) || (fileid < 3)) {
  } else {
    cst = fclose(filestar);
    if (cst == 0) {
      f = 0;
      fileid = (signed char)(fileid - 2);
      eml_openfiles[fileid - 1] = NULL;
      eml_autoflush[fileid - 1] = true;
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char e_fileManager()
{
  signed char f;
  signed char j;
  char cv4[9];
  int i1;
  static const char cv5[9] = { 'N', 'o', 'd', 'e', '.', 'i', 'm', 'g', '\x00' };

  char cv6[3];
  static const char cv7[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i1 = 0; i1 < 9; i1++) {
      cv4[i1] = cv5[i1];
    }

    for (i1 = 0; i1 < 3; i1++) {
      cv6[i1] = cv7[i1];
    }

    filestar = fopen(cv4, cv6);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : double varargin_1
// Return Type  : FILE *
//
FILE * f_fileManager(double varargin_1)
{
  FILE * f;
  signed char fileid;
  fileid = (signed char)rt_roundd_snf(varargin_1);
  if ((fileid > 22) || (fileid < 0) || (varargin_1 != fileid)) {
    fileid = -1;
  }

  if (fileid >= 3) {
    f = eml_openfiles[fileid - 3];
  } else if (fileid == 0) {
    f = stdin;
  } else if (fileid == 1) {
    f = stdout;
  } else if (fileid == 2) {
    f = stderr;
  } else {
    f = NULL;
  }

  return f;
}

//
// Arguments    : FILE * *f
//                boolean_T *a
// Return Type  : void
//
void fileManager(FILE * *f, boolean_T *a)
{
  *f = stdout;
  *a = true;
}

//
// Arguments    : void
// Return Type  : void
//
void filedata_init()
{
  FILE * a;
  int i;
  a = NULL;
  for (i = 0; i < 20; i++) {
    eml_autoflush[i] = false;
    eml_openfiles[i] = a;
  }
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char g_fileManager()
{
  signed char f;
  signed char j;
  char cv8[9];
  int i2;
  static const char cv9[9] = { 'E', 'l', 'e', 'm', '.', 'i', 'm', 'g', '\x00' };

  char cv10[3];
  static const char cv11[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i2 = 0; i2 < 9; i2++) {
      cv8[i2] = cv9[i2];
    }

    for (i2 = 0; i2 < 3; i2++) {
      cv10[i2] = cv11[i2];
    }

    filestar = fopen(cv8, cv10);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char h_fileManager()
{
  signed char f;
  signed char j;
  char cv12[9];
  int i3;
  static const char cv13[9] = { 'F', 'a', 'c', 'e', '.', 'i', 'm', 'g', '\x00' };

  char cv14[3];
  static const char cv15[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i3 = 0; i3 < 9; i3++) {
      cv12[i3] = cv13[i3];
    }

    for (i3 = 0; i3 < 3; i3++) {
      cv14[i3] = cv15[i3];
    }

    filestar = fopen(cv12, cv14);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : const emxArray_char_T *varargin_1
// Return Type  : signed char
//
signed char i_fileManager(const emxArray_char_T *varargin_1)
{
  signed char f;
  signed char j;
  emxArray_char_T *r4;
  int i20;
  int loop_ub;
  char cv20[3];
  static const char cv21[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    emxInit_char_T(&r4, 2);
    i20 = r4->size[0] * r4->size[1];
    r4->size[0] = 1;
    r4->size[1] = varargin_1->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)r4, i20, (int)sizeof(char));
    loop_ub = varargin_1->size[1];
    for (i20 = 0; i20 < loop_ub; i20++) {
      r4->data[r4->size[0] * i20] = varargin_1->data[varargin_1->size[0] * i20];
    }

    r4->data[r4->size[0] * varargin_1->size[1]] = '\x00';
    for (i20 = 0; i20 < 3; i20++) {
      cv20[i20] = cv21[i20];
    }

    filestar = fopen(&r4->data[0], cv20);
    emxFree_char_T(&r4);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : const emxArray_char_T *varargin_1
// Return Type  : signed char
//
signed char j_fileManager(const emxArray_char_T *varargin_1)
{
  signed char f;
  signed char j;
  emxArray_char_T *r7;
  int i24;
  int loop_ub;
  char cv25[3];
  static const char cv26[3] = { 'w', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    emxInit_char_T(&r7, 2);
    i24 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = varargin_1->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)r7, i24, (int)sizeof(char));
    loop_ub = varargin_1->size[1];
    for (i24 = 0; i24 < loop_ub; i24++) {
      r7->data[r7->size[0] * i24] = varargin_1->data[varargin_1->size[0] * i24];
    }

    r7->data[r7->size[0] * varargin_1->size[1]] = '\x00';
    for (i24 = 0; i24 < 3; i24++) {
      cv25[i24] = cv26[i24];
    }

    filestar = fopen(&r7->data[0], cv25);
    emxFree_char_T(&r7);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char k_fileManager()
{
  signed char f;
  signed char j;
  char cv28[37];
  int i26;
  static const char cv29[37] = { 'm', 'n', 'i', '_', 'i', 'c', 'b', 'm', '1',
    '5', '2', '_', 't', '1', '_', 't', 'a', 'l', '_', 'n', 'l', 'i', 'n', '_',
    'a', 's', 'y', 'm', '_', '0', '9', 'c', '.', 'h', 'd', 'r', '\x00' };

  char cv30[3];
  static const char cv31[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i26 = 0; i26 < 37; i26++) {
      cv28[i26] = cv29[i26];
    }

    for (i26 = 0; i26 < 3; i26++) {
      cv30[i26] = cv31[i26];
    }

    filestar = fopen(cv28, cv30);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char l_fileManager()
{
  signed char f;
  signed char j;
  char cv32[37];
  int i27;
  static const char cv33[37] = { 'm', 'n', 'i', '_', 'i', 'c', 'b', 'm', '1',
    '5', '2', '_', 't', '1', '_', 't', 'a', 'l', '_', 'n', 'l', 'i', 'n', '_',
    'a', 's', 'y', 'm', '_', '0', '9', 'c', '.', 'i', 'm', 'g', '\x00' };

  char cv34[3];
  static const char cv35[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i27 = 0; i27 < 37; i27++) {
      cv32[i27] = cv33[i27];
    }

    for (i27 = 0; i27 < 3; i27++) {
      cv34[i27] = cv35[i27];
    }

    filestar = fopen(cv32, cv34);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char m_fileManager()
{
  signed char f;
  signed char j;
  char cv36[9];
  int i28;
  static const char cv37[9] = { 'g', 'r', 'a', 'y', '.', 'i', 'm', 'g', '\x00' };

  char cv38[3];
  static const char cv39[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i28 = 0; i28 < 9; i28++) {
      cv36[i28] = cv37[i28];
    }

    for (i28 = 0; i28 < 3; i28++) {
      cv38[i28] = cv39[i28];
    }

    filestar = fopen(cv36, cv38);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char n_fileManager()
{
  signed char f;
  signed char j;
  char cv40[10];
  int i29;
  static const char cv41[10] = { 'w', 'h', 'i', 't', 'e', '.', 'i', 'm', 'g',
    '\x00' };

  char cv42[3];
  static const char cv43[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i29 = 0; i29 < 10; i29++) {
      cv40[i29] = cv41[i29];
    }

    for (i29 = 0; i29 < 3; i29++) {
      cv42[i29] = cv43[i29];
    }

    filestar = fopen(cv40, cv42);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

//
// Arguments    : void
// Return Type  : signed char
//
signed char o_fileManager()
{
  signed char f;
  signed char j;
  char cv44[8];
  int i30;
  static const char cv45[8] = { 'c', 's', 'f', '.', 'i', 'm', 'g', '\x00' };

  char cv46[3];
  static const char cv47[3] = { 'r', 'b', '\x00' };

  FILE * filestar;
  f = -1;
  j = filedata();
  if (j < 1) {
  } else {
    for (i30 = 0; i30 < 8; i30++) {
      cv44[i30] = cv45[i30];
    }

    for (i30 = 0; i30 < 3; i30++) {
      cv46[i30] = cv47[i30];
    }

    filestar = fopen(cv44, cv46);
    if (filestar != (FILE *)NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      f = (signed char)(j + 2);
    }
  }

  return f;
}

