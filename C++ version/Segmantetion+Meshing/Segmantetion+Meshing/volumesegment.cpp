
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
#include "imdilate.h"
#include "fclose.h"
#include "fread.h"
#include "fopen.h"
#include "prod.h"
#include <stdio.h>


// Function Declarations
static void c_eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y);

// Function Definitions

//
// Arguments    : const emxArray_boolean_T *x
//                emxArray_int32_T *y
// Return Type  : void
//
static void c_eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y)
{
  int n;
  int k;
  int i;
  int j;
  n = x->size[0] * x->size[1] * x->size[2];
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      k++;
    }
  }

  j = y->size[0];
  y->size[0] = k;
  emxEnsureCapacity((emxArray__common *)y, j, (int)sizeof(int));
  j = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      y->data[j] = i;
      j++;
    }
  }
}

//
// volumesegment segments an anatomical MRI.
// Arguments    : const double mri_dim[3]
//                const emxArray_real_T *mri_anatomy
//                emxArray_real_T *segmented
// Return Type  : void
//
void b_volumesegment(const double mri_dim[3], const emxArray_real_T *mri_anatomy,
                     emxArray_real_T *segmented)
{
  emxArray_real_T *V1;
  double fid;
  int iv5[3];
  int braindil;
  emxArray_real_T *b_V1;
  int k;
  double img_siz;
  emxArray_real_T *V2;
  emxArray_real_T *V3;
  emxArray_boolean_T *scalpmask;
  emxArray_boolean_T *brainmask;
  emxArray_boolean_T *b_brainmask;
  emxArray_boolean_T *b_braindil;
  int c_braindil;
  emxArray_boolean_T *c_brainmask;
  emxArray_int32_T *r11;
  emxArray_boolean_T *d_braindil;
  int b_scalpmask;
  int c_scalpmask;
  int d_scalpmask;
  int e_scalpmask;
  emxArray_boolean_T *d_brainmask;
  emxArray_boolean_T *cleanbrain;
  emxArray_boolean_T *e_braindil;
  emxInit_real_T(&V1, 1);

  
  fid = k_fopen();
  r_fread(fid, prod(mri_dim), V1);

  // 'int16'
  for (braindil = 0; braindil < 3; braindil++) {
    iv5[braindil] = (int)mri_dim[braindil];
  }

  c_emxInit_real_T(&b_V1, 3);
  braindil = b_V1->size[0] * b_V1->size[1] * b_V1->size[2];
  b_V1->size[0] = iv5[0];
  b_V1->size[1] = iv5[1];
  b_V1->size[2] = iv5[2];
  emxEnsureCapacity((emxArray__common *)b_V1, braindil, (int)sizeof(double));
  for (k = 0; k + 1 <= V1->size[0]; k++) {
    b_V1->data[k] = V1->data[k];
  }

  b_fclose(fid);
  fid = l_fopen();
  img_siz = mri_dim[0];
  for (k = 0; k < 2; k++) {
    img_siz *= mri_dim[k + 1];
  }

  r_fread(fid, img_siz, V1);

  // 'int16'
  for (braindil = 0; braindil < 3; braindil++) {
    iv5[braindil] = (int)mri_dim[braindil];
  }

  c_emxInit_real_T(&V2, 3);
  braindil = V2->size[0] * V2->size[1] * V2->size[2];
  V2->size[0] = iv5[0];
  V2->size[1] = iv5[1];
  V2->size[2] = iv5[2];
  emxEnsureCapacity((emxArray__common *)V2, braindil, (int)sizeof(double));
  for (k = 0; k + 1 <= V1->size[0]; k++) {
    V2->data[k] = V1->data[k];
  }

  b_fclose(fid);
  fid = m_fopen();
  img_siz = mri_dim[0];
  for (k = 0; k < 2; k++) {
    img_siz *= mri_dim[k + 1];
  }

  r_fread(fid, img_siz, V1);

  // 'int16'
  for (braindil = 0; braindil < 3; braindil++) {
    iv5[braindil] = (int)mri_dim[braindil];
  }

  c_emxInit_real_T(&V3, 3);
  braindil = V3->size[0] * V3->size[1] * V3->size[2];
  V3->size[0] = iv5[0];
  V3->size[1] = iv5[1];
  V3->size[2] = iv5[2];
  emxEnsureCapacity((emxArray__common *)V3, braindil, (int)sizeof(double));
  for (k = 0; k + 1 <= V1->size[0]; k++) {
    V3->data[k] = V1->data[k];
  }

  emxFree_real_T(&V1);
  b_emxInit_boolean_T(&scalpmask, 7);
  c_emxInit_boolean_T(&brainmask, 3);
  b_fclose(fid);

  //  create the requested output fields
  //  create scalpmask - no tpm or brainmask is required to create it
  e_fprintf();

  //  fill the slices along each dimension (because using a single one is
  //  just arbitrary, and behavior depends on how the voxeldata is in the
  //  volume.
  b_fillholes(mri_anatomy, scalpmask);

  //  threshold again to remove little parts outside of head
  //  create the brain from the tpm
  g_fprintf();

  //  make binary mask from brain
  braindil = brainmask->size[0] * brainmask->size[1] * brainmask->size[2];
  brainmask->size[0] = b_V1->size[0];
  brainmask->size[1] = b_V1->size[1];
  brainmask->size[2] = b_V1->size[2];
  emxEnsureCapacity((emxArray__common *)brainmask, braindil, (int)sizeof
                    (boolean_T));
  k = b_V1->size[0] * b_V1->size[1] * b_V1->size[2];
  for (braindil = 0; braindil < k; braindil++) {
    brainmask->data[braindil] = ((b_V1->data[braindil] + V2->data[braindil]) +
      V3->data[braindil] > 0.0);
  }

  emxFree_real_T(&V3);
  emxFree_real_T(&V2);
  emxFree_real_T(&b_V1);
  c_emxInit_boolean_T(&b_brainmask, 3);

  //  output: brain
  //  create skull from brain mask FIXME check this (e.g. strel_bol)
  i_fprintf();
  braindil = b_brainmask->size[0] * b_brainmask->size[1] * b_brainmask->size[2];
  b_brainmask->size[0] = brainmask->size[0];
  b_brainmask->size[1] = brainmask->size[1];
  b_brainmask->size[2] = brainmask->size[2];
  emxEnsureCapacity((emxArray__common *)b_brainmask, braindil, (int)sizeof
                    (boolean_T));
  k = brainmask->size[0] * brainmask->size[1] * brainmask->size[2];
  for (braindil = 0; braindil < k; braindil++) {
    b_brainmask->data[braindil] = ((int)brainmask->data[braindil] > 0);
  }

  c_emxInit_boolean_T(&b_braindil, 3);
  imdilate(b_brainmask, b_braindil);
  braindil = b_braindil->size[0] * b_braindil->size[1] * b_braindil->size[2];
  emxEnsureCapacity((emxArray__common *)b_braindil, braindil, (int)sizeof
                    (boolean_T));
  k = b_braindil->size[0];
  braindil = b_braindil->size[1];
  c_braindil = b_braindil->size[2];
  k = k * braindil * c_braindil;
  emxFree_boolean_T(&b_brainmask);
  for (braindil = 0; braindil < k; braindil++) {
    b_braindil->data[braindil] = (b_braindil->data[braindil] &&
      (!brainmask->data[braindil]));
  }

  c_emxInit_boolean_T(&c_brainmask, 3);
  braindil = c_brainmask->size[0] * c_brainmask->size[1] * c_brainmask->size[2];
  c_brainmask->size[0] = brainmask->size[0];
  c_brainmask->size[1] = brainmask->size[1];
  c_brainmask->size[2] = brainmask->size[2];
  emxEnsureCapacity((emxArray__common *)c_brainmask, braindil, (int)sizeof
                    (boolean_T));
  k = brainmask->size[0] * brainmask->size[1] * brainmask->size[2];
  for (braindil = 0; braindil < k; braindil++) {
    c_brainmask->data[braindil] = ((int)brainmask->data[braindil] > 0);
  }

  b_emxInit_int32_T(&r11, 1);
  c_eml_li_find(c_brainmask, r11);
  k = r11->size[0];
  emxFree_boolean_T(&c_brainmask);
  for (braindil = 0; braindil < k; braindil++) {
    scalpmask->data[r11->data[braindil] - 1] = false;
  }

  c_emxInit_boolean_T(&d_braindil, 3);
  braindil = d_braindil->size[0] * d_braindil->size[1] * d_braindil->size[2];
  d_braindil->size[0] = b_braindil->size[0];
  d_braindil->size[1] = b_braindil->size[1];
  d_braindil->size[2] = b_braindil->size[2];
  emxEnsureCapacity((emxArray__common *)d_braindil, braindil, (int)sizeof
                    (boolean_T));
  k = b_braindil->size[0] * b_braindil->size[1] * b_braindil->size[2];
  for (braindil = 0; braindil < k; braindil++) {
    d_braindil->data[braindil] = ((int)b_braindil->data[braindil] > 0);
  }

  c_eml_li_find(d_braindil, r11);
  k = r11->size[0];
  emxFree_boolean_T(&d_braindil);
  for (braindil = 0; braindil < k; braindil++) {
    scalpmask->data[r11->data[braindil] - 1] = false;
  }

  emxFree_int32_T(&r11);
  k_fprintf();

  //  fill holes in the head image and create the canonical binary volume
  k_fprintf();
  braindil = scalpmask->size[0] * scalpmask->size[1] * scalpmask->size[2] *
    scalpmask->size[3] * scalpmask->size[4] * scalpmask->size[5] *
    scalpmask->size[6];
  emxEnsureCapacity((emxArray__common *)scalpmask, braindil, (int)sizeof
                    (boolean_T));
  k = scalpmask->size[0];
  braindil = scalpmask->size[1];
  c_braindil = scalpmask->size[2];
  b_scalpmask = scalpmask->size[3];
  c_scalpmask = scalpmask->size[4];
  d_scalpmask = scalpmask->size[5];
  e_scalpmask = scalpmask->size[6];
  k = k * braindil * c_braindil * b_scalpmask * c_scalpmask * d_scalpmask *
    e_scalpmask;
  for (braindil = 0; braindil < k; braindil++) {
    scalpmask->data[braindil] = ((int)scalpmask->data[braindil] > 0);
  }

  c_emxInit_boolean_T(&d_brainmask, 3);
  braindil = d_brainmask->size[0] * d_brainmask->size[1] * d_brainmask->size[2];
  d_brainmask->size[0] = brainmask->size[0];
  d_brainmask->size[1] = brainmask->size[1];
  d_brainmask->size[2] = brainmask->size[2];
  emxEnsureCapacity((emxArray__common *)d_brainmask, braindil, (int)sizeof
                    (boolean_T));
  k = brainmask->size[0] * brainmask->size[1] * brainmask->size[2];
  for (braindil = 0; braindil < k; braindil++) {
    d_brainmask->data[braindil] = ((int)brainmask->data[braindil] > 0);
  }

  c_emxInit_boolean_T(&cleanbrain, 3);
  c_emxInit_boolean_T(&e_braindil, 3);
  c_fillholes(d_brainmask, cleanbrain);
  braindil = e_braindil->size[0] * e_braindil->size[1] * e_braindil->size[2];
  e_braindil->size[0] = b_braindil->size[0];
  e_braindil->size[1] = b_braindil->size[1];
  e_braindil->size[2] = b_braindil->size[2];
  emxEnsureCapacity((emxArray__common *)e_braindil, braindil, (int)sizeof
                    (boolean_T));
  k = b_braindil->size[0] * b_braindil->size[1] * b_braindil->size[2];
  emxFree_boolean_T(&d_brainmask);
  for (braindil = 0; braindil < k; braindil++) {
    e_braindil->data[braindil] = ((int)b_braindil->data[braindil] > 0);
  }

  emxFree_boolean_T(&b_braindil);
  c_fillholes(e_braindil, brainmask);

  //  add brain image as additional segment
  braindil = cleanbrain->size[0] * cleanbrain->size[1] * cleanbrain->size[2];
  emxEnsureCapacity((emxArray__common *)cleanbrain, braindil, (int)sizeof
                    (boolean_T));
  k = cleanbrain->size[0];
  braindil = cleanbrain->size[1];
  c_braindil = cleanbrain->size[2];
  k = k * braindil * c_braindil;
  emxFree_boolean_T(&e_braindil);
  for (braindil = 0; braindil < k; braindil++) {
    cleanbrain->data[braindil] = ((int)cleanbrain->data[braindil] > 0);
  }

  k = scalpmask->size[0];
  braindil = segmented->size[0] * segmented->size[1] * segmented->size[2];
  segmented->size[0] = k;
  emxEnsureCapacity((emxArray__common *)segmented, braindil, (int)sizeof(double));
  k = scalpmask->size[1];
  braindil = segmented->size[0] * segmented->size[1] * segmented->size[2];
  segmented->size[1] = k;
  emxEnsureCapacity((emxArray__common *)segmented, braindil, (int)sizeof(double));
  k = scalpmask->size[2];
  braindil = segmented->size[0] * segmented->size[1] * segmented->size[2];
  segmented->size[2] = k;
  emxEnsureCapacity((emxArray__common *)segmented, braindil, (int)sizeof(double));
  k = scalpmask->size[0] * scalpmask->size[1] * scalpmask->size[2];
  for (braindil = 0; braindil < k; braindil++) {
    segmented->data[braindil] = ((double)scalpmask->data[braindil] + (double)
      cleanbrain->data[braindil]) + (double)((int)brainmask->data[braindil] > 0);
  }

  emxFree_boolean_T(&cleanbrain);
  emxFree_boolean_T(&brainmask);
  emxFree_boolean_T(&scalpmask);

  
}

//
// volumesegment segments an anatomical MRI.
// Arguments    : const struct1_T *mri
//                double segmented[8530021]
// Return Type  : void
//
void volumesegment(const struct1_T *mri, double segmented[8530021])
{
  emxArray_real_T *V1;
  double fid;
  int iv6[3];
  int k;
  emxArray_real_T *b_V1;
  double img_siz;
  emxArray_real_T *V2;
  emxArray_real_T *V3;
  emxArray_boolean_T *brainmask;
  static boolean_T scalpmask[8530021];
  int loop_ub;
  emxArray_boolean_T *b_brainmask;
  emxArray_boolean_T *braindil;
  int b_braindil;
  emxArray_boolean_T *c_brainmask;
  emxArray_int32_T *r15;
  emxArray_boolean_T *c_braindil;
  emxArray_boolean_T *d_brainmask;
  emxArray_boolean_T *cleanbrain;
  emxArray_boolean_T *d_braindil;
  emxInit_real_T(&V1, 1);

  
  fid = k_fopen();
  r_fread(fid, prod(mri->dim), V1);

  // 'int16'
  for (k = 0; k < 3; k++) {
    iv6[k] = (int)mri->dim[k];
  }

  c_emxInit_real_T(&b_V1, 3);
  k = b_V1->size[0] * b_V1->size[1] * b_V1->size[2];
  b_V1->size[0] = iv6[0];
  b_V1->size[1] = iv6[1];
  b_V1->size[2] = iv6[2];
  emxEnsureCapacity((emxArray__common *)b_V1, k, (int)sizeof(double));
  for (k = 0; k + 1 <= V1->size[0]; k++) {
    b_V1->data[k] = V1->data[k];
  }

  b_fclose(fid);
  fid = l_fopen();
  img_siz = mri->dim[0];
  for (k = 0; k < 2; k++) {
    img_siz *= mri->dim[k + 1];
  }

  r_fread(fid, img_siz, V1);

  // 'int16'
  for (k = 0; k < 3; k++) {
    iv6[k] = (int)mri->dim[k];
  }

  c_emxInit_real_T(&V2, 3);
  k = V2->size[0] * V2->size[1] * V2->size[2];
  V2->size[0] = iv6[0];
  V2->size[1] = iv6[1];
  V2->size[2] = iv6[2];
  emxEnsureCapacity((emxArray__common *)V2, k, (int)sizeof(double));
  for (k = 0; k + 1 <= V1->size[0]; k++) {
    V2->data[k] = V1->data[k];
  }

  b_fclose(fid);
  fid = m_fopen();
  img_siz = mri->dim[0];
  for (k = 0; k < 2; k++) {
    img_siz *= mri->dim[k + 1];
  }

  r_fread(fid, img_siz, V1);

  // 'int16'
  for (k = 0; k < 3; k++) {
    iv6[k] = (int)mri->dim[k];
  }

  c_emxInit_real_T(&V3, 3);
  k = V3->size[0] * V3->size[1] * V3->size[2];
  V3->size[0] = iv6[0];
  V3->size[1] = iv6[1];
  V3->size[2] = iv6[2];
  emxEnsureCapacity((emxArray__common *)V3, k, (int)sizeof(double));
  for (k = 0; k + 1 <= V1->size[0]; k++) {
    V3->data[k] = V1->data[k];
  }

  emxFree_real_T(&V1);
  c_emxInit_boolean_T(&brainmask, 3);
  b_fclose(fid);

  //  create the requested output fields
  //  create scalpmask - no tpm or brainmask is required to create it
  e_fprintf();

  //  fill the slices along each dimension (because using a single one is
  //  just arbitrary, and behavior depends on how the voxeldata is in the
  //  volume.

  fillholes(mri->anatomy, scalpmask);

  //  threshold again to remove little parts outside of head
  //  create the brain from the tpm
  g_fprintf();

  //  make binary mask from brain
  k = brainmask->size[0] * brainmask->size[1] * brainmask->size[2];
  brainmask->size[0] = b_V1->size[0];
  brainmask->size[1] = b_V1->size[1];
  brainmask->size[2] = b_V1->size[2];
  emxEnsureCapacity((emxArray__common *)brainmask, k, (int)sizeof(boolean_T));
  loop_ub = b_V1->size[0] * b_V1->size[1] * b_V1->size[2];
  for (k = 0; k < loop_ub; k++) {
    brainmask->data[k] = ((b_V1->data[k] + V2->data[k]) + V3->data[k] > 0.0);
  }

  emxFree_real_T(&V3);
  emxFree_real_T(&V2);
  emxFree_real_T(&b_V1);
  c_emxInit_boolean_T(&b_brainmask, 3);

  //  output: brain
  //  create skull from brain mask FIXME check this (e.g. strel_bol)
  i_fprintf();
  k = b_brainmask->size[0] * b_brainmask->size[1] * b_brainmask->size[2];
  b_brainmask->size[0] = brainmask->size[0];
  b_brainmask->size[1] = brainmask->size[1];
  b_brainmask->size[2] = brainmask->size[2];
  emxEnsureCapacity((emxArray__common *)b_brainmask, k, (int)sizeof(boolean_T));
  loop_ub = brainmask->size[0] * brainmask->size[1] * brainmask->size[2];
  for (k = 0; k < loop_ub; k++) {
    b_brainmask->data[k] = ((int)brainmask->data[k] > 0);
  }

  c_emxInit_boolean_T(&braindil, 3);
  imdilate(b_brainmask, braindil);
  k = braindil->size[0] * braindil->size[1] * braindil->size[2];
  emxEnsureCapacity((emxArray__common *)braindil, k, (int)sizeof(boolean_T));
  k = braindil->size[0];
  loop_ub = braindil->size[1];
  b_braindil = braindil->size[2];
  loop_ub = k * loop_ub * b_braindil;
  emxFree_boolean_T(&b_brainmask);
  for (k = 0; k < loop_ub; k++) {
    braindil->data[k] = (braindil->data[k] && (!brainmask->data[k]));
  }

  c_emxInit_boolean_T(&c_brainmask, 3);
  k = c_brainmask->size[0] * c_brainmask->size[1] * c_brainmask->size[2];
  c_brainmask->size[0] = brainmask->size[0];
  c_brainmask->size[1] = brainmask->size[1];
  c_brainmask->size[2] = brainmask->size[2];
  emxEnsureCapacity((emxArray__common *)c_brainmask, k, (int)sizeof(boolean_T));
  loop_ub = brainmask->size[0] * brainmask->size[1] * brainmask->size[2];
  for (k = 0; k < loop_ub; k++) {
    c_brainmask->data[k] = ((int)brainmask->data[k] > 0);
  }

  b_emxInit_int32_T(&r15, 1);
  c_eml_li_find(c_brainmask, r15);
  loop_ub = r15->size[0];
  emxFree_boolean_T(&c_brainmask);
  for (k = 0; k < loop_ub; k++) {
    scalpmask[r15->data[k] - 1] = false;
  }

  c_emxInit_boolean_T(&c_braindil, 3);
  k = c_braindil->size[0] * c_braindil->size[1] * c_braindil->size[2];
  c_braindil->size[0] = braindil->size[0];
  c_braindil->size[1] = braindil->size[1];
  c_braindil->size[2] = braindil->size[2];
  emxEnsureCapacity((emxArray__common *)c_braindil, k, (int)sizeof(boolean_T));
  loop_ub = braindil->size[0] * braindil->size[1] * braindil->size[2];
  for (k = 0; k < loop_ub; k++) {
    c_braindil->data[k] = ((int)braindil->data[k] > 0);
  }

  c_eml_li_find(c_braindil, r15);
  loop_ub = r15->size[0];
  emxFree_boolean_T(&c_braindil);
  for (k = 0; k < loop_ub; k++) {
    scalpmask[r15->data[k] - 1] = false;
  }

  emxFree_int32_T(&r15);
  c_emxInit_boolean_T(&d_brainmask, 3);
  k_fprintf();

  //  fill holes in the head image and create the canonical binary volume
  k_fprintf();
  k = d_brainmask->size[0] * d_brainmask->size[1] * d_brainmask->size[2];
  d_brainmask->size[0] = brainmask->size[0];
  d_brainmask->size[1] = brainmask->size[1];
  d_brainmask->size[2] = brainmask->size[2];
  emxEnsureCapacity((emxArray__common *)d_brainmask, k, (int)sizeof(boolean_T));
  loop_ub = brainmask->size[0] * brainmask->size[1] * brainmask->size[2];
  for (k = 0; k < loop_ub; k++) {
    d_brainmask->data[k] = ((int)brainmask->data[k] > 0);
  }

  c_emxInit_boolean_T(&cleanbrain, 3);
  c_emxInit_boolean_T(&d_braindil, 3);
  c_fillholes(d_brainmask, cleanbrain);
  k = d_braindil->size[0] * d_braindil->size[1] * d_braindil->size[2];
  d_braindil->size[0] = braindil->size[0];
  d_braindil->size[1] = braindil->size[1];
  d_braindil->size[2] = braindil->size[2];
  emxEnsureCapacity((emxArray__common *)d_braindil, k, (int)sizeof(boolean_T));
  loop_ub = braindil->size[0] * braindil->size[1] * braindil->size[2];
  emxFree_boolean_T(&d_brainmask);
  for (k = 0; k < loop_ub; k++) {
    d_braindil->data[k] = ((int)braindil->data[k] > 0);
  }

  emxFree_boolean_T(&braindil);
  c_fillholes(d_braindil, brainmask);

  //  add brain image as additional segment
  emxFree_boolean_T(&d_braindil);
  for (k = 0; k < 8530021; k++) {
    segmented[k] = ((double)((int)scalpmask[k] > 0) + (double)((int)
      cleanbrain->data[k] > 0)) + (double)((int)brainmask->data[k] > 0);
  }

  emxFree_boolean_T(&cleanbrain);
  emxFree_boolean_T(&brainmask);

 
}
