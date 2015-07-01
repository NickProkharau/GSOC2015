
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
#include "fprintf.h"
#include "fileManager.h"
#include <stdio.h>


// Function Declarations
static double f_fprintf();
static double h_fprintf();
static double j_fprintf();

// Function Definitions

//
// Arguments    : void
// Return Type  : double
//
static double f_fprintf()
{
  int nbytesint;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[20] = { 'c', 'r', 'e', 'a', 't', 'i', 'n', 'g', ' ',
    's', 'c', 'a', 'l', 'p', 'm', 'a', 's', 'k', '\x0a', '\x00' };

  nbytesint = 0;
  fileManager(&filestar, &autoflush);
  if (filestar == NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt);
    if (autoflush) {
      fflush(filestar);
    }
  }

  return nbytesint;
}

//
// Arguments    : void
// Return Type  : double
//
static double h_fprintf()
{
  int nbytesint;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[20] = { 'c', 'r', 'e', 'a', 't', 'i', 'n', 'g', ' ',
    'b', 'r', 'a', 'i', 'n', 'm', 'a', 's', 'k', '\x0a', '\x00' };

  nbytesint = 0;
  fileManager(&filestar, &autoflush);
  if (filestar == NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt);
    if (autoflush) {
      fflush(filestar);
    }
  }

  return nbytesint;
}

//
// Arguments    : void
// Return Type  : double
//
static double j_fprintf()
{
  int nbytesint;
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[20] = { 'c', 'r', 'e', 'a', 't', 'i', 'n', 'g', ' ',
    's', 'k', 'u', 'l', 'l', 'm', 'a', 's', 'k', '\x0a', '\x00' };

  nbytesint = 0;
  fileManager(&filestar, &autoflush);
  if (filestar == NULL) {
  } else {
    nbytesint = fprintf(filestar, cfmt);
    if (autoflush) {
      fflush(filestar);
    }
  }

  return nbytesint;
}

//
// Arguments    : void
// Return Type  : void
//
void b_fprintf()
{
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[70] = { 'c', 'r', 'e', 'a', 't', 'i', 'n', 'g', ' ',
    's', 'u', 'r', 'f', 'a', 'c', 'e', ' ', 'a', 'n', 'd', ' ', 't', 'e', 't',
    'r', 'a', 'h', 'e', 'd', 'r', 'a', 'l', ' ', 'm', 'e', 's', 'h', ' ', 'f',
    'r', 'o', 'm', ' ', 'a', ' ', 'm', 'u', 'l', 't', 'i', '-', 'd', 'o', 'm',
    'a', 'i', 'n', ' ', 'v', 'o', 'l', 'u', 'm', 'e', ' ', '.', '.', '.', '\x0a',
    '\x00' };

  fileManager(&filestar, &autoflush);
  if (filestar == NULL) {
  } else {
    fprintf(filestar, cfmt);
    if (autoflush) {
      fflush(filestar);
    }
  }
}

//
// Arguments    : double varargin_4
// Return Type  : void
//
void c_fprintf(double varargin_4)
{
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[58] = { 'n', 'o', 'd', 'e', ' ', 'n', 'u', 'm', 'b',
    'e', 'r', ':', '	', '%', 'f', '\x0a', 't', 'r', 'i', 'a', 'n', 'g', 'l', 'e',
    's', ':', '	', '%', 'f', '\x0a', 't', 'e', 't', 'r', 'a', 'h', 'e', 'd', 'r',
    'a', ':', '	', '%', 'f', '\x0a', 'r', 'e', 'g', 'i', 'o', 'n', 's', ':', '	',
    '%', 'f', '\x0a', '\x00' };

  fileManager(&filestar, &autoflush);
  if (filestar == NULL) {
  } else {
    fprintf(filestar, cfmt, 25141.0, 90592.0, 113150.0, varargin_4);
    if (autoflush) {
      fflush(filestar);
    }
  }
}

//
// Arguments    : void
// Return Type  : void
//
void d_fprintf()
{
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[36] = { 's', 'u', 'r', 'f', 'a', 'c', 'e', ' ', 'a',
    'n', 'd', ' ', 'v', 'o', 'l', 'u', 'm', 'e', ' ', 'm', 'e', 's', 'h', 'e',
    's', ' ', 'c', 'o', 'm', 'p', 'l', 'e', 't', 'e', '\x0a', '\x00' };

  fileManager(&filestar, &autoflush);
  if (filestar == NULL) {
  } else {
    fprintf(filestar, cfmt);
    if (autoflush) {
      fflush(filestar);
    }
  }
}

//
// Arguments    : void
// Return Type  : void
//
void e_fprintf()
{
  f_fprintf();
}

//
// Arguments    : void
// Return Type  : void
//
void g_fprintf()
{
  h_fprintf();
}

//
// Arguments    : void
// Return Type  : void
//
void i_fprintf()
{
  j_fprintf();
}

//
// Arguments    : void
// Return Type  : void
//
void k_fprintf()
{
  boolean_T autoflush;
  FILE * filestar;
  static const char cfmt[43] = { 'f', 'i', 'l', 'l', 'i', 'n', 'g', ' ', 'h',
    'o', 'l', 'e', 's', ' ', 'i', 'n', ' ', 't', 'h', 'e', ' ', 'v', 'o', 'l',
    'u', 'm', 'e', 't', 'r', 'i', 'c', ' ', 'i', 'm', 'a', 'g', 'e', 's', '.',
    '.', '.', '\x0a', '\x00' };

  fileManager(&filestar, &autoflush);
  if (filestar == NULL) {
  } else {
    fprintf(filestar, cfmt);
    if (autoflush) {
      fflush(filestar);
    }
  }
}

