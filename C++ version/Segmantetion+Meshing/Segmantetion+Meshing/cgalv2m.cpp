
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
#include "Segmentation1_emxutil.h"
#include "fprintf.h"
#include "fclose.h"
#include "fopen.h"
#include <stdio.h>


// Function Definitions

//
// convert a binary (or multi-valued) volume to tetrahedral mesh
// Arguments    : const unsigned char vol[8530021]
//                double opt
//                double maxvol
//                double node[100564]
//                double elem[565750]
//                double face[362368]
// Return Type  : void
//
void cgalv2m(const unsigned char [8530021], double, double, double node[100564],
             double elem[565750], double face[362368])
{
  double fid;
  static int idx[113150];
  int k;
  boolean_T p;
  static int idx0[113150];
  int i;
  int i2;
  int j;
  int pEnd;
  int nb;
  int b_k;
  int qEnd;
  int kEnd;
  emxArray_real_T *b;
  double x;
  int32_T exitg1;
  int exponent;
  emxArray_real_T *b_b;


  b_fprintf();

  //
  //  saveinr(vol,fname)
  //
  //  save a surface mesh to INR Format
  //
  //  input:
  //       vol: input, a binary volume
  //       fname: output file name
  fid = b_fopen();

  // fwrite(fid,header,'char');
  // fwrite(fid,vol,dtype);
  b_fclose(fid);

  //  "U+623F U+9A9E"
  // pwd
  // system(cmd);
  b_readmedit(node, elem, face);

  // node=node*double(vol(1,1,1));
  for (k = 0; k < 113150; k++) {
    idx[k] = k + 1;
  }

  for (k = 0; k < 113150; k += 2) {
    if ((elem[452600 + k] <= elem[k + 452601]) || rtIsNaN(elem[k + 452601])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  for (i = 0; i < 113150; i++) {
    idx0[i] = 1;
  }

  i = 2;
  while (i < 113150) {
    i2 = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 113151; pEnd = qEnd + i) {
      nb = j;
      b_k = pEnd - 1;
      qEnd = j + i2;
      if (qEnd > 113151) {
        qEnd = 113151;
      }

      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        if ((elem[idx[nb - 1] + 452599] <= elem[idx[b_k] + 452599]) || rtIsNaN
            (elem[idx[b_k] + 452599])) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          idx0[k] = idx[nb - 1];
          nb++;
          if (nb == pEnd) {
            while (b_k + 1 < qEnd) {
              k++;
              idx0[k] = idx[b_k];
              b_k++;
            }
          }
        } else {
          idx0[k] = idx[b_k];
          b_k++;
          if (b_k + 1 == qEnd) {
            while (nb < pEnd) {
              k++;
              idx0[k] = idx[nb - 1];
              nb++;
            }
          }
        }

        k++;
      }

      for (k = 0; k + 1 <= kEnd; k++) {
        idx[(j + k) - 1] = idx0[k];
      }

      j = qEnd;
    }

    i = i2;
  }

  emxInit_real_T(&b, 1);
  i2 = b->size[0];
  b->size[0] = 113150;
  emxEnsureCapacity((emxArray__common *)b, i2, (int)sizeof(double));
  for (k = 0; k < 113150; k++) {
    b->data[k] = elem[idx[k] + 452599];
  }

  k = 0;
  while ((k + 1 <= 113150) && rtIsInf(b->data[k]) && (b->data[k] < 0.0)) {
    k++;
  }

  b_k = k;
  k = 113150;
  while ((k >= 1) && rtIsNaN(b->data[k - 1])) {
    k--;
  }

  pEnd = 113150 - k;
  while ((k >= 1) && rtIsInf(b->data[k - 1]) && (b->data[k - 1] > 0.0)) {
    k--;
  }

  i2 = 113150 - (k + pEnd);
  nb = -1;
  if (b_k > 0) {
    nb = 0;
  }

  i = (b_k + k) - b_k;
  while (b_k + 1 <= i) {
    x = b->data[b_k];
    do {
      exitg1 = 0;
      b_k++;
      if (b_k + 1 > i) {
        exitg1 = 1;
      } else {
        fid = fabs(x / 2.0);
        if ((!rtIsInf(fid)) && (!rtIsNaN(fid))) {
          if (fid <= 2.2250738585072014E-308) {
            fid = 4.94065645841247E-324;
          } else {
            frexp(fid, &exponent);
            fid = ldexp(1.0, exponent - 53);
          }
        } else {
          fid = rtNaN;
        }

        if ((fabs(x - b->data[b_k]) < fid) || (rtIsInf(b->data[b_k]) && rtIsInf
             (x) && ((b->data[b_k] > 0.0) == (x > 0.0)))) {
          p = true;
        } else {
          p = false;
        }

        if (!p) {
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);

    nb++;
    b->data[nb] = x;
  }

  if (i2 > 0) {
    nb++;
    b->data[nb] = b->data[i];
  }

  b_k = i + i2;
  for (j = 1; j <= pEnd; j++) {
    nb++;
    b->data[nb] = b->data[(b_k + j) - 1];
  }

  if (1 > nb + 1) {
    i = -1;
  } else {
    i = nb;
  }

  emxInit_real_T(&b_b, 1);
  i2 = b_b->size[0];
  b_b->size[0] = i + 1;
  emxEnsureCapacity((emxArray__common *)b_b, i2, (int)sizeof(double));
  for (i2 = 0; i2 <= i; i2++) {
    b_b->data[i2] = b->data[i2];
  }

  i2 = b->size[0];
  b->size[0] = b_b->size[0];
  emxEnsureCapacity((emxArray__common *)b, i2, (int)sizeof(double));
  i = b_b->size[0];
  for (i2 = 0; i2 < i; i2++) {
    b->data[i2] = b_b->data[i2];
  }

  emxFree_real_T(&b_b);
  c_fprintf((double)b->size[0]);
  d_fprintf();
  emxFree_real_T(&b);
}


