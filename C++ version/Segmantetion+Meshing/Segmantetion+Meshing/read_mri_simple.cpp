
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
#include "norm.h"
#include "Segmentation1_emxutil.h"
#include "eml_sort.h"
#include "fread.h"
#include "fopen.h"
#include "diag.h"
#include "eye.h"
#include "fclose.h"
#include "Segmentation1_rtwutil.h"
#include <stdio.h>


// Function Declarations
static void b_eml_li_find(const boolean_T x[8], emxArray_int32_T *y);
static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y);
static void estimate_units(double size, char unit[2]);
static void idrange(const double x[24], double r[3]);

// Function Definitions

//
// Arguments    : const boolean_T x[8]
//                emxArray_int32_T *y
// Return Type  : void
//
static void b_eml_li_find(const boolean_T x[8], emxArray_int32_T *y)
{
  int k;
  int i;
  int j;
  k = 0;
  for (i = 0; i < 8; i++) {
    if (x[i]) {
      k++;
    }
  }

  j = y->size[0];
  y->size[0] = k;
  emxEnsureCapacity((emxArray__common *)y, j, (int)sizeof(int));
  j = 0;
  for (i = 0; i < 8; i++) {
    if (x[i]) {
      y->data[j] = i + 1;
      j++;
    }
  }
}

//
// Arguments    : const emxArray_boolean_T *x
//                emxArray_int32_T *y
// Return Type  : void
//
static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y)
{
  int k;
  int i;
  int j;
  k = 0;
  for (i = 1; i <= x->size[1]; i++) {
    if (x->data[i - 1]) {
      k++;
    }
  }

  j = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = k;
  emxEnsureCapacity((emxArray__common *)y, j, (int)sizeof(int));
  j = 0;
  for (i = 1; i <= x->size[1]; i++) {
    if (x->data[i - 1]) {
      y->data[j] = i;
      j++;
    }
  }
}

//
// Arguments    : double size
//                char unit[2]
// Return Type  : void
//
static void estimate_units(double size, char unit[2])
{
  double indx;
  int i23;
  static const char cv22[2] = { 'm', ' ' };

  static const char cv23[2] = { 'd', 'm' };

  static const char cv24[2] = { 'c', 'm' };

  indx = rt_roundd_snf((double)(log10(size) + 1.8));
  if (indx > 4.0) {
    indx = 4.0;

    // warning('assuming that the units are "%s"', 'mm');
  } else {
    if (indx < 1.0) {
      indx = 1.0;

      // warning('assuming that the units are "%s"', 'm');
    }
  }

  switch ((int)indx) {
   case 1:
    for (i23 = 0; i23 < 2; i23++) {
      unit[i23] = cv22[i23];
    }
    break;

   case 2:
    for (i23 = 0; i23 < 2; i23++) {
      unit[i23] = cv23[i23];
    }
    break;

   case 3:
    for (i23 = 0; i23 < 2; i23++) {
      unit[i23] = cv24[i23];
    }
    break;

   case 4:
    for (i23 = 0; i23 < 2; i23++) {
      unit[i23] = 'm';
    }
    break;

   default:
    for (i23 = 0; i23 < 2; i23++) {
      unit[i23] = cv22[i23];
    }
    break;
  }
}

//
// Arguments    : const double x[24]
//                double r[3]
// Return Type  : void
//
static void idrange(const double x[24], double r[3])
{
  int i;
  boolean_T keeprow[8];
  int b_r;
  emxArray_real_T *sx;
  emxArray_int32_T *r5;
  int iyStart;
  emxArray_real_T *b_sx;
  emxArray_real_T *r6;
  unsigned int varargin_2_idx_1;
  unsigned int ii[2];
  double c_r;
  double d0;
  for (i = 0; i < 8; i++) {
    keeprow[i] = true;
  }

  for (i = 0; i < 3; i++) {
    for (b_r = 0; b_r < 8; b_r++) {
      keeprow[b_r] = (keeprow[b_r] && ((!rtIsInf(x[b_r + (i << 3)])) &&
        (!rtIsNaN(x[b_r + (i << 3)]))));
    }
  }

  b_emxInit_real_T(&sx, 2);
  b_emxInit_int32_T(&r5, 1);
  b_eml_li_find(keeprow, r5);
  b_r = sx->size[0] * sx->size[1];
  sx->size[0] = r5->size[0];
  sx->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)sx, b_r, (int)sizeof(double));
  for (b_r = 0; b_r < 3; b_r++) {
    i = r5->size[0];
    for (iyStart = 0; iyStart < i; iyStart++) {
      sx->data[iyStart + sx->size[0] * b_r] = x[(r5->data[iyStart] + (b_r << 3))
        - 1];
    }
  }

  b_emxInit_real_T(&b_sx, 2);
  b_r = b_sx->size[0] * b_sx->size[1];
  b_sx->size[0] = sx->size[0];
  b_sx->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)b_sx, b_r, (int)sizeof(double));
  i = sx->size[0] * sx->size[1];
  for (b_r = 0; b_r < i; b_r++) {
    b_sx->data[b_r] = sx->data[b_r];
  }

  b_emxInit_real_T(&r6, 2);
  eml_sort(b_sx, sx);
  b_eml_li_find(keeprow, r5);
  b_r = r6->size[0] * r6->size[1];
  r6->size[0] = r5->size[0];
  r6->size[1] = 3;
  emxEnsureCapacity((emxArray__common *)r6, b_r, (int)sizeof(double));
  emxFree_real_T(&b_sx);
  for (b_r = 0; b_r < 3; b_r++) {
    i = r5->size[0];
    for (iyStart = 0; iyStart < i; iyStart++) {
      r6->data[iyStart + r6->size[0] * b_r] = x[(r5->data[iyStart] + (b_r << 3))
        - 1];
    }
  }

  emxFree_int32_T(&r5);
  varargin_2_idx_1 = (unsigned int)r6->size[0];
  emxFree_real_T(&r6);
  for (b_r = 0; b_r < 2; b_r++) {
    c_r = 0.1 + 0.8 * (double)b_r;
    if (1U == varargin_2_idx_1) {
      d0 = 1.0;
    } else {
      d0 = (1.0 - c_r) + c_r * (double)varargin_2_idx_1;
    }

    ii[b_r] = (unsigned int)rt_roundd_snf(d0);
  }

  //  indices for 10 & 90 percentile
  i = 1;
  iyStart = 0;
  for (b_r = 0; b_r < 3; b_r++) {
    r[iyStart] = sx->data[((int)ii[i % 2] + sx->size[0] * (i / 2)) - 1] -
      sx->data[((int)ii[(i - 1) % 2] + sx->size[0] * ((i - 1) / 2)) - 1];
    i += 2;
    iyStart++;
  }

  emxFree_real_T(&sx);
}

//
// read_mri_simple reads anatomical and functional MRI data from different
//  file formats. The output data is structured in such a way that it is
//  comparable
// Arguments    : double mri_dim[3]
//                emxArray_real_T *mri_anatomy
//                double mri_transform[16]
//                char mri_unit[2]
//                char mri_coordsys[8]
// Return Type  : void
//
void b_read_mri_simple(double mri_dim[3], emxArray_real_T *mri_anatomy, double
  mri_transform[16], char mri_unit[2], char mri_coordsys[8])
{
  emxArray_int16_T *t0_dim;
  emxArray_int16_T *r8;
  emxArray_int32_T *fp;
  emxArray_uint8_T *b_fp;
  emxArray_uint8_T *c_fp;
  emxArray_int32_T *d_fp;
  emxArray_int16_T *e_fp;
  emxArray_uint8_T *f_fp;
  emxArray_uint8_T *g_fp;
  double h_fp;
  int i25;
  int n;
  emxArray_real32_T *t0_pixdim;
  emxArray_real32_T *r9;
  emxArray_real32_T *i_fp;
  emxArray_real32_T *j_fp;
  emxArray_real32_T *k_fp;
  emxArray_int16_T *l_fp;
  emxArray_int16_T *m_fp;
  emxArray_int16_T *n_fp;
  emxArray_int16_T *o_fp;
  emxArray_real32_T *t0_quatern_b;
  emxArray_real32_T *p_fp;
  emxArray_real32_T *q_fp;
  emxArray_real32_T *r_fp;
  emxArray_int16_T *s_fp;
  emxArray_uint8_T *t_fp;
  emxArray_uint8_T *u_fp;
  emxArray_real32_T *v_fp;
  emxArray_real32_T *w_fp;
  emxArray_real32_T *x_fp;
  emxArray_real32_T *y_fp;
  emxArray_int32_T *ab_fp;
  emxArray_int32_T *bb_fp;
  emxArray_uint8_T *cb_fp;
  emxArray_uint8_T *db_fp;
  emxArray_int16_T *eb_fp;
  emxArray_int16_T *fb_fp;
  emxArray_real32_T *t0_quatern_c;
  emxArray_real32_T *t0_quatern_d;
  emxArray_real32_T *t0_qoffset_x;
  emxArray_real32_T *t0_qoffset_y;
  emxArray_real32_T *t0_qoffset_z;
  emxArray_real_T *dim;
  emxArray_real32_T *gb_fp;
  emxArray_real32_T *hb_fp;
  emxArray_real32_T *ib_fp;
  emxArray_uint8_T *jb_fp;
  emxArray_uint8_T *kb_fp;
  emxArray_real_T *pixdim;
  emxArray_real_T *b_t0_quatern_b;
  emxArray_real_T *img;
  double R[16];
  emxArray_real_T *b_img;
  double dv3[12];
  int ar;
  double T[16];
  emxArray_real_T *Z;
  int unnamed_idx_1;
  int br;
  emxArray_boolean_T *b_Z;
  emxArray_int32_T *r10;
  emxArray_real_T *c_Z;
  double y[16];
  emxArray_real_T *M;
  unsigned int varargin_1_idx_1;
  int ic;
  int ib;
  int ia;
  static const signed char b[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1, -1,
    -1, 1 };

  static const signed char iv4[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1,
    -1, -1, 1 };

  double img_siz;
  short varargin_1[7];
  unsigned short uv1[7];
  emxArray_real_T *c_img;
  static const char cv27[8] = { 'n', 'e', 'u', 'r', 'o', 'm', 'a', 'g' };

  double pos_head[24];
  double dv4[3];
  emxInit_int16_T(&t0_dim, 2);
  b_emxInit_int16_T(&r8, 1);
  b_emxInit_int32_T(&fp, 1);
  emxInit_uint8_T(&b_fp, 1);
  emxInit_uint8_T(&c_fp, 1);
  b_emxInit_int32_T(&d_fp, 1);
  b_emxInit_int16_T(&e_fp, 1);
  emxInit_uint8_T(&f_fp, 1);
  emxInit_uint8_T(&g_fp, 1);

  //  Open it if possible
  h_fp = i_fopen();

  //  Read the fields
  e_fread(h_fp, fp);
  f_fread(h_fp, b_fp);
  g_fread(h_fp, c_fp);
  e_fread(h_fp, d_fp);
  h_fread(h_fp, e_fp);
  i_fread(h_fp, f_fp);
  i_fread(h_fp, g_fp);
  j_fread(h_fp, r8);
  i25 = t0_dim->size[0] * t0_dim->size[1];
  t0_dim->size[0] = 1;
  t0_dim->size[1] = r8->size[0];
  emxEnsureCapacity((emxArray__common *)t0_dim, i25, (int)sizeof(short));
  n = r8->size[0];
  emxFree_uint8_T(&g_fp);
  emxFree_uint8_T(&f_fp);
  emxFree_int16_T(&e_fp);
  emxFree_int32_T(&d_fp);
  emxFree_uint8_T(&c_fp);
  emxFree_uint8_T(&b_fp);
  emxFree_int32_T(&fp);
  for (i25 = 0; i25 < n; i25++) {
    t0_dim->data[t0_dim->size[0] * i25] = r8->data[i25];
  }

  emxFree_int16_T(&r8);
  emxInit_real32_T(&t0_pixdim, 2);
  b_emxInit_real32_T(&r9, 1);
  b_emxInit_real32_T(&i_fp, 1);
  b_emxInit_real32_T(&j_fp, 1);
  b_emxInit_real32_T(&k_fp, 1);
  b_emxInit_int16_T(&l_fp, 1);
  b_emxInit_int16_T(&m_fp, 1);
  b_emxInit_int16_T(&n_fp, 1);
  b_emxInit_int16_T(&o_fp, 1);
  k_fread(h_fp, i_fp);
  k_fread(h_fp, j_fp);
  k_fread(h_fp, k_fp);
  h_fread(h_fp, l_fp);
  h_fread(h_fp, m_fp);
  h_fread(h_fp, n_fp);
  h_fread(h_fp, o_fp);
  l_fread(h_fp, r9);
  i25 = t0_pixdim->size[0] * t0_pixdim->size[1];
  t0_pixdim->size[0] = 1;
  t0_pixdim->size[1] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)t0_pixdim, i25, (int)sizeof(float));
  n = r9->size[0];
  emxFree_int16_T(&o_fp);
  emxFree_int16_T(&n_fp);
  emxFree_int16_T(&m_fp);
  emxFree_int16_T(&l_fp);
  emxFree_real32_T(&k_fp);
  emxFree_real32_T(&j_fp);
  emxFree_real32_T(&i_fp);
  for (i25 = 0; i25 < n; i25++) {
    t0_pixdim->data[t0_pixdim->size[0] * i25] = r9->data[i25];
  }

  emxInit_real32_T(&t0_quatern_b, 2);
  b_emxInit_real32_T(&p_fp, 1);
  b_emxInit_real32_T(&q_fp, 1);
  b_emxInit_real32_T(&r_fp, 1);
  b_emxInit_int16_T(&s_fp, 1);
  emxInit_uint8_T(&t_fp, 1);
  emxInit_uint8_T(&u_fp, 1);
  b_emxInit_real32_T(&v_fp, 1);
  b_emxInit_real32_T(&w_fp, 1);
  b_emxInit_real32_T(&x_fp, 1);
  b_emxInit_real32_T(&y_fp, 1);
  b_emxInit_int32_T(&ab_fp, 1);
  b_emxInit_int32_T(&bb_fp, 1);
  emxInit_uint8_T(&cb_fp, 1);
  emxInit_uint8_T(&db_fp, 1);
  b_emxInit_int16_T(&eb_fp, 1);
  b_emxInit_int16_T(&fb_fp, 1);
  k_fread(h_fp, p_fp);
  k_fread(h_fp, q_fp);
  k_fread(h_fp, r_fp);
  h_fread(h_fp, s_fp);
  i_fread(h_fp, t_fp);
  i_fread(h_fp, u_fp);
  k_fread(h_fp, v_fp);
  k_fread(h_fp, w_fp);
  k_fread(h_fp, x_fp);
  k_fread(h_fp, y_fp);
  e_fread(h_fp, ab_fp);
  e_fread(h_fp, bb_fp);
  m_fread(h_fp, cb_fp);
  n_fread(h_fp, db_fp);
  h_fread(h_fp, eb_fp);
  h_fread(h_fp, fb_fp);
  k_fread(h_fp, r9);
  i25 = t0_quatern_b->size[0] * t0_quatern_b->size[1];
  t0_quatern_b->size[0] = 1;
  t0_quatern_b->size[1] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)t0_quatern_b, i25, (int)sizeof(float));
  n = r9->size[0];
  emxFree_int16_T(&fb_fp);
  emxFree_int16_T(&eb_fp);
  emxFree_uint8_T(&db_fp);
  emxFree_uint8_T(&cb_fp);
  emxFree_int32_T(&bb_fp);
  emxFree_int32_T(&ab_fp);
  emxFree_real32_T(&y_fp);
  emxFree_real32_T(&x_fp);
  emxFree_real32_T(&w_fp);
  emxFree_real32_T(&v_fp);
  emxFree_uint8_T(&u_fp);
  emxFree_uint8_T(&t_fp);
  emxFree_int16_T(&s_fp);
  emxFree_real32_T(&r_fp);
  emxFree_real32_T(&q_fp);
  emxFree_real32_T(&p_fp);
  for (i25 = 0; i25 < n; i25++) {
    t0_quatern_b->data[t0_quatern_b->size[0] * i25] = r9->data[i25];
  }

  emxInit_real32_T(&t0_quatern_c, 2);
  k_fread(h_fp, r9);
  i25 = t0_quatern_c->size[0] * t0_quatern_c->size[1];
  t0_quatern_c->size[0] = 1;
  t0_quatern_c->size[1] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)t0_quatern_c, i25, (int)sizeof(float));
  n = r9->size[0];
  for (i25 = 0; i25 < n; i25++) {
    t0_quatern_c->data[t0_quatern_c->size[0] * i25] = r9->data[i25];
  }

  emxInit_real32_T(&t0_quatern_d, 2);
  k_fread(h_fp, r9);
  i25 = t0_quatern_d->size[0] * t0_quatern_d->size[1];
  t0_quatern_d->size[0] = 1;
  t0_quatern_d->size[1] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)t0_quatern_d, i25, (int)sizeof(float));
  n = r9->size[0];
  for (i25 = 0; i25 < n; i25++) {
    t0_quatern_d->data[t0_quatern_d->size[0] * i25] = r9->data[i25];
  }

  emxInit_real32_T(&t0_qoffset_x, 2);
  k_fread(h_fp, r9);
  i25 = t0_qoffset_x->size[0] * t0_qoffset_x->size[1];
  t0_qoffset_x->size[0] = 1;
  t0_qoffset_x->size[1] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)t0_qoffset_x, i25, (int)sizeof(float));
  n = r9->size[0];
  for (i25 = 0; i25 < n; i25++) {
    t0_qoffset_x->data[t0_qoffset_x->size[0] * i25] = r9->data[i25];
  }

  emxInit_real32_T(&t0_qoffset_y, 2);
  k_fread(h_fp, r9);
  i25 = t0_qoffset_y->size[0] * t0_qoffset_y->size[1];
  t0_qoffset_y->size[0] = 1;
  t0_qoffset_y->size[1] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)t0_qoffset_y, i25, (int)sizeof(float));
  n = r9->size[0];
  for (i25 = 0; i25 < n; i25++) {
    t0_qoffset_y->data[t0_qoffset_y->size[0] * i25] = r9->data[i25];
  }

  emxInit_real32_T(&t0_qoffset_z, 2);
  k_fread(h_fp, r9);
  i25 = t0_qoffset_z->size[0] * t0_qoffset_z->size[1];
  t0_qoffset_z->size[0] = 1;
  t0_qoffset_z->size[1] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)t0_qoffset_z, i25, (int)sizeof(float));
  n = r9->size[0];
  for (i25 = 0; i25 < n; i25++) {
    t0_qoffset_z->data[t0_qoffset_z->size[0] * i25] = r9->data[i25];
  }

  emxFree_real32_T(&r9);
  b_emxInit_real_T(&dim, 2);
  b_emxInit_real32_T(&gb_fp, 1);
  b_emxInit_real32_T(&hb_fp, 1);
  b_emxInit_real32_T(&ib_fp, 1);
  emxInit_uint8_T(&jb_fp, 1);
  emxInit_uint8_T(&kb_fp, 1);
  o_fread(h_fp, gb_fp);
  o_fread(h_fp, hb_fp);
  o_fread(h_fp, ib_fp);
  p_fread(h_fp, jb_fp);
  q_fread(h_fp, kb_fp);
  b_fclose(h_fp);
  i25 = dim->size[0] * dim->size[1];
  dim->size[0] = 1;
  dim->size[1] = t0_dim->size[1];
  emxEnsureCapacity((emxArray__common *)dim, i25, (int)sizeof(double));
  n = t0_dim->size[0] * t0_dim->size[1];
  emxFree_uint8_T(&kb_fp);
  emxFree_uint8_T(&jb_fp);
  emxFree_real32_T(&ib_fp);
  emxFree_real32_T(&hb_fp);
  emxFree_real32_T(&gb_fp);
  for (i25 = 0; i25 < n; i25++) {
    dim->data[i25] = t0_dim->data[i25];
  }

  emxFree_int16_T(&t0_dim);
  b_emxInit_real_T(&pixdim, 2);

  // Decode qform info from NIFTI-1 headers
  i25 = pixdim->size[0] * pixdim->size[1];
  pixdim->size[0] = 1;
  pixdim->size[1] = t0_pixdim->size[1];
  emxEnsureCapacity((emxArray__common *)pixdim, i25, (int)sizeof(double));
  n = t0_pixdim->size[0] * t0_pixdim->size[1];
  for (i25 = 0; i25 < n; i25++) {
    pixdim->data[i25] = t0_pixdim->data[i25];
  }

  emxFree_real32_T(&t0_pixdim);
  b_emxInit_real_T(&b_t0_quatern_b, 2);

  //  Rotations from quaternions
  i25 = b_t0_quatern_b->size[0] * b_t0_quatern_b->size[1];
  b_t0_quatern_b->size[0] = 1;
  b_t0_quatern_b->size[1] = (t0_quatern_b->size[1] + t0_quatern_c->size[1]) +
    t0_quatern_d->size[1];
  emxEnsureCapacity((emxArray__common *)b_t0_quatern_b, i25, (int)sizeof(double));
  n = t0_quatern_b->size[1];
  for (i25 = 0; i25 < n; i25++) {
    b_t0_quatern_b->data[b_t0_quatern_b->size[0] * i25] = t0_quatern_b->
      data[t0_quatern_b->size[0] * i25];
  }

  n = t0_quatern_c->size[1];
  for (i25 = 0; i25 < n; i25++) {
    b_t0_quatern_b->data[b_t0_quatern_b->size[0] * (i25 + t0_quatern_b->size[1])]
      = t0_quatern_c->data[t0_quatern_c->size[0] * i25];
  }

  n = t0_quatern_d->size[1];
  for (i25 = 0; i25 < n; i25++) {
    b_t0_quatern_b->data[b_t0_quatern_b->size[0] * ((i25 + t0_quatern_b->size[1])
      + t0_quatern_c->size[1])] = t0_quatern_d->data[t0_quatern_d->size[0] * i25];
  }

  emxFree_real32_T(&t0_quatern_d);
  emxFree_real32_T(&t0_quatern_c);
  emxFree_real32_T(&t0_quatern_b);
  emxInit_real_T(&img, 1);
  b_rot_matrix(b_t0_quatern_b, R);

  //  Translations
  i25 = img->size[0];
  img->size[0] = ((t0_qoffset_x->size[1] + t0_qoffset_y->size[1]) +
                  t0_qoffset_z->size[1]) + 1;
  emxEnsureCapacity((emxArray__common *)img, i25, (int)sizeof(double));
  n = t0_qoffset_x->size[1];
  emxFree_real_T(&b_t0_quatern_b);
  for (i25 = 0; i25 < n; i25++) {
    img->data[i25] = t0_qoffset_x->data[t0_qoffset_x->size[0] * i25];
  }

  n = t0_qoffset_y->size[1];
  for (i25 = 0; i25 < n; i25++) {
    img->data[i25 + t0_qoffset_x->size[1]] = t0_qoffset_y->data
      [t0_qoffset_y->size[0] * i25];
  }

  n = t0_qoffset_z->size[1];
  for (i25 = 0; i25 < n; i25++) {
    img->data[(i25 + t0_qoffset_x->size[1]) + t0_qoffset_y->size[1]] =
      t0_qoffset_z->data[t0_qoffset_z->size[0] * i25];
  }

  b_emxInit_real_T(&b_img, 2);
  img->data[(t0_qoffset_x->size[1] + t0_qoffset_y->size[1]) + t0_qoffset_z->
    size[1]] = 1.0;
  eye(dv3);
  n = img->size[0];
  i25 = b_img->size[0] * b_img->size[1];
  b_img->size[0] = n;
  b_img->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)b_img, i25, (int)sizeof(double));
  emxFree_real32_T(&t0_qoffset_z);
  emxFree_real32_T(&t0_qoffset_y);
  emxFree_real32_T(&t0_qoffset_x);
  for (i25 = 0; i25 < n; i25++) {
    b_img->data[i25] = img->data[i25];
  }

  for (i25 = 0; i25 < 3; i25++) {
    for (ar = 0; ar < 4; ar++) {
      T[ar + (i25 << 2)] = dv3[ar + (i25 << 2)];
    }
  }

  for (i25 = 0; i25 < 4; i25++) {
    T[12 + i25] = b_img->data[i25];
  }

  emxFree_real_T(&b_img);

  //  Zooms.  Note that flips are derived from the first
  //  element of pixdim, which is normally unused.
  if (dim->data[0] <= 3.0) {
    n = (int)dim->data[0] + 1;
  } else {
    n = 4;
  }

  if (2 > n) {
    i25 = 1;
    ar = 1;
  } else {
    i25 = 2;
    ar = n + 1;
  }

  b_emxInit_real_T(&Z, 2);
  unnamed_idx_1 = 5 - n;
  br = Z->size[0] * Z->size[1];
  Z->size[0] = 1;
  Z->size[1] = (ar - i25) + unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)Z, br, (int)sizeof(double));
  n = ar - i25;
  for (br = 0; br < n; br++) {
    Z->data[Z->size[0] * br] = pixdim->data[(i25 + br) - 1];
  }

  for (br = 0; br < unnamed_idx_1; br++) {
    Z->data[Z->size[0] * ((br + ar) - i25)] = 1.0;
  }

  emxInit_boolean_T(&b_Z, 2);
  i25 = b_Z->size[0] * b_Z->size[1];
  b_Z->size[0] = 1;
  b_Z->size[1] = Z->size[1];
  emxEnsureCapacity((emxArray__common *)b_Z, i25, (int)sizeof(boolean_T));
  n = Z->size[0] * Z->size[1];
  for (i25 = 0; i25 < n; i25++) {
    b_Z->data[i25] = (Z->data[i25] < 0.0);
  }

  emxInit_int32_T(&r10, 2);
  eml_li_find(b_Z, r10);
  n = r10->size[0] * r10->size[1];
  emxFree_boolean_T(&b_Z);
  for (i25 = 0; i25 < n; i25++) {
    Z->data[r10->data[i25] - 1] = 1.0;
  }

  emxFree_int32_T(&r10);
  if (pixdim->data[0] < 0.0) {
    Z->data[2] = -Z->data[2];
  }

  emxFree_real_T(&pixdim);
  b_emxInit_real_T(&c_Z, 2);
  diag(Z, c_Z);
  emxFree_real_T(&Z);
  for (i25 = 0; i25 < 4; i25++) {
    for (ar = 0; ar < 4; ar++) {
      y[i25 + (ar << 2)] = 0.0;
      for (br = 0; br < 4; br++) {
        y[i25 + (ar << 2)] += T[i25 + (br << 2)] * R[br + (ar << 2)];
      }
    }
  }

  b_emxInit_real_T(&M, 2);
  if (c_Z->size[0] == 1) {
    i25 = M->size[0] * M->size[1];
    M->size[0] = 4;
    M->size[1] = c_Z->size[1];
    emxEnsureCapacity((emxArray__common *)M, i25, (int)sizeof(double));
    for (i25 = 0; i25 < 4; i25++) {
      n = c_Z->size[1];
      for (ar = 0; ar < n; ar++) {
        M->data[i25 + M->size[0] * ar] = 0.0;
        for (br = 0; br < 4; br++) {
          M->data[i25 + M->size[0] * ar] += y[i25 + (br << 2)] * c_Z->data[br +
            c_Z->size[0] * ar];
        }
      }
    }
  } else {
    varargin_1_idx_1 = (unsigned int)c_Z->size[1];
    i25 = M->size[0] * M->size[1];
    M->size[0] = 4;
    emxEnsureCapacity((emxArray__common *)M, i25, (int)sizeof(double));
    i25 = M->size[0] * M->size[1];
    M->size[1] = (int)varargin_1_idx_1;
    emxEnsureCapacity((emxArray__common *)M, i25, (int)sizeof(double));
    n = (int)varargin_1_idx_1 << 2;
    for (i25 = 0; i25 < n; i25++) {
      M->data[i25] = 0.0;
    }

    if (c_Z->size[1] == 0) {
    } else {
      n = (c_Z->size[1] - 1) << 2;
      for (unnamed_idx_1 = 0; unnamed_idx_1 <= n; unnamed_idx_1 += 4) {
        for (ic = unnamed_idx_1 + 1; ic <= unnamed_idx_1 + 4; ic++) {
          M->data[ic - 1] = 0.0;
        }
      }

      br = 0;
      for (unnamed_idx_1 = 0; unnamed_idx_1 <= n; unnamed_idx_1 += 4) {
        ar = 0;
        for (ib = br; ib + 1 <= br + 4; ib++) {
          if (c_Z->data[ib] != 0.0) {
            ia = ar;
            for (ic = unnamed_idx_1; ic + 1 <= unnamed_idx_1 + 4; ic++) {
              ia++;
              M->data[ic] += c_Z->data[ib] * y[ia - 1];
            }
          }

          ar += 4;
        }

        br += 4;
      }
    }
  }

  emxFree_real_T(&c_Z);

  //  Convert from first voxel at [1,1,1]
  //  to first voxel at [0,0,0]
  if (M->size[1] == 1) {
    for (i25 = 0; i25 < 4; i25++) {
      for (ar = 0; ar < 4; ar++) {
        R[i25 + (ar << 2)] = 0.0;
        for (br = 0; br < 4; br++) {
          R[i25 + (ar << 2)] += M->data[i25 + (br << 2)] * (double)b[br + (ar <<
            2)];
        }
      }
    }
  } else {
    memset(&R[0], 0, sizeof(double) << 4);
    for (unnamed_idx_1 = 0; unnamed_idx_1 < 14; unnamed_idx_1 += 4) {
      for (ic = unnamed_idx_1; ic + 1 <= unnamed_idx_1 + 4; ic++) {
        R[ic] = 0.0;
      }
    }

    br = 0;
    for (unnamed_idx_1 = 0; unnamed_idx_1 < 14; unnamed_idx_1 += 4) {
      ar = 0;
      i25 = br + M->size[1];
      for (ib = br; ib + 1 <= i25; ib++) {
        if (iv4[ib] != 0) {
          ia = ar;
          for (ic = unnamed_idx_1; ic + 1 <= unnamed_idx_1 + 4; ic++) {
            ia++;
            R[ic] += (double)iv4[ib] * M->data[ia - 1];
          }
        }

        ar += 4;
      }

      br += M->size[1];
    }
  }

  emxFree_real_T(&M);

  // read volumes
  h_fp = j_fopen();

  // fseek(fid, 0, 'bof');
  img_siz = dim->data[1];
  for (n = 0; n < 6; n++) {
    img_siz *= dim->data[n + 2];
  }

  r_fread(h_fp, img_siz, img);

  // 'int16'
  for (i25 = 0; i25 < 7; i25++) {
    varargin_1[i25] = (short)dim->data[i25 + 1];
  }

  emxFree_real_T(&dim);
  for (i25 = 0; i25 < 7; i25++) {
    uv1[i25] = (unsigned short)varargin_1[i25];
  }

  d_emxInit_real_T(&c_img, 7);
  i25 = c_img->size[0] * c_img->size[1] * c_img->size[2] * c_img->size[3] *
    c_img->size[4] * c_img->size[5] * c_img->size[6];
  c_img->size[0] = uv1[0];
  c_img->size[1] = uv1[1];
  c_img->size[2] = uv1[2];
  c_img->size[3] = uv1[3];
  c_img->size[4] = uv1[4];
  c_img->size[5] = uv1[5];
  c_img->size[6] = uv1[6];
  emxEnsureCapacity((emxArray__common *)c_img, i25, (int)sizeof(double));
  for (n = 0; n + 1 <= img->size[0]; n++) {
    c_img->data[n] = img->data[n];
  }

  emxFree_real_T(&img);

  //  set up the axes of the volume in voxel coordinates
  mri_dim[0] = (unsigned short)c_img->size[0];
  mri_dim[1] = (unsigned short)c_img->size[1];
  mri_dim[2] = (unsigned short)c_img->size[2];
  i25 = mri_anatomy->size[0] * mri_anatomy->size[1] * mri_anatomy->size[2] *
    mri_anatomy->size[3] * mri_anatomy->size[4] * mri_anatomy->size[5] *
    mri_anatomy->size[6];
  mri_anatomy->size[0] = c_img->size[0];
  mri_anatomy->size[1] = c_img->size[1];
  mri_anatomy->size[2] = c_img->size[2];
  mri_anatomy->size[3] = c_img->size[3];
  mri_anatomy->size[4] = c_img->size[4];
  mri_anatomy->size[5] = c_img->size[5];
  mri_anatomy->size[6] = c_img->size[6];
  emxEnsureCapacity((emxArray__common *)mri_anatomy, i25, (int)sizeof(double));
  n = c_img->size[0] * c_img->size[1] * c_img->size[2] * c_img->size[3] *
    c_img->size[4] * c_img->size[5] * c_img->size[6];
  for (i25 = 0; i25 < n; i25++) {
    mri_anatomy->data[i25] = c_img->data[i25];
  }

  emxFree_real_T(&c_img);
  memcpy(&mri_transform[0], &R[0], sizeof(double) << 4);
  for (i25 = 0; i25 < 8; i25++) {
    mri_coordsys[i25] = cv27[i25];
  }

  //  try to determine the units of the coordinate system
  corner_points(mri_dim, R, pos_head);
  idrange(pos_head, dv4);
  estimate_units(norm(dv4), mri_unit);
}

//
// read_mri_simple reads anatomical and functional MRI data from different
//  file formats. The output data is structured in such a way that it is
//  comparable
// Arguments    : const emxArray_char_T *fname
//                struct0_T *mri
// Return Type  : void
//
void read_mri_simple(const emxArray_char_T *fname, struct0_T *mri)
{
  emxArray_char_T *iname;
  int i19;
  int in_idx;
  int copyfrom;
  int out_idx;
  int b_in_idx;
  int b_copyfrom;
  int out_len;
  int patt_idx;
  int tmp_in_idx;
  static const char cv16[4] = { '.', 'h', 'd', 'r' };

  static const char cv17[4] = { '.', 'i', 'm', 'g' };

  emxArray_int16_T *t1_dim;
  emxArray_int16_T *r1;
  emxArray_int32_T *fp;
  emxArray_uint8_T *b_fp;
  emxArray_uint8_T *c_fp;
  emxArray_int32_T *d_fp;
  emxArray_int16_T *e_fp;
  emxArray_uint8_T *f_fp;
  emxArray_uint8_T *g_fp;
  double h_fp;
  emxArray_real32_T *t1_pixdim;
  emxArray_real32_T *r2;
  emxArray_real32_T *i_fp;
  emxArray_real32_T *j_fp;
  emxArray_real32_T *k_fp;
  emxArray_int16_T *l_fp;
  emxArray_int16_T *m_fp;
  emxArray_int16_T *n_fp;
  emxArray_int16_T *o_fp;
  emxArray_real32_T *t1_quatern_b;
  emxArray_real32_T *p_fp;
  emxArray_real32_T *q_fp;
  emxArray_real32_T *r_fp;
  emxArray_int16_T *s_fp;
  emxArray_uint8_T *t_fp;
  emxArray_uint8_T *u_fp;
  emxArray_real32_T *v_fp;
  emxArray_real32_T *w_fp;
  emxArray_real32_T *x_fp;
  emxArray_real32_T *y_fp;
  emxArray_int32_T *ab_fp;
  emxArray_int32_T *bb_fp;
  emxArray_uint8_T *cb_fp;
  emxArray_uint8_T *db_fp;
  emxArray_int16_T *eb_fp;
  emxArray_int16_T *fb_fp;
  emxArray_real32_T *t1_quatern_c;
  emxArray_real32_T *t1_quatern_d;
  emxArray_real32_T *t1_qoffset_x;
  emxArray_real32_T *t1_qoffset_y;
  emxArray_real32_T *t1_qoffset_z;
  emxArray_real_T *dim;
  emxArray_real32_T *gb_fp;
  emxArray_real32_T *hb_fp;
  emxArray_real32_T *ib_fp;
  emxArray_uint8_T *jb_fp;
  emxArray_uint8_T *kb_fp;
  emxArray_real_T *pixdim;
  emxArray_real_T *b_t1_quatern_b;
  emxArray_real_T *img;
  double R[16];
  emxArray_real_T *b_img;
  double dv1[12];
  double T[16];
  emxArray_real_T *Z;
  emxArray_boolean_T *b_Z;
  emxArray_int32_T *r3;
  emxArray_real_T *c_Z;
  double y[16];
  emxArray_real_T *M;
  unsigned int varargin_1_idx_1;
  static const signed char b[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1, -1,
    -1, 1 };

  static const signed char iv1[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1,
    -1, -1, 1 };

  double img_siz;
  short varargin_1[7];
  unsigned short uv0[7];
  emxArray_real_T *c_img;
  static const char cv18[8] = { 'n', 'e', 'u', 'r', 'o', 'm', 'a', 'g' };

  double pos_head[24];
  double dv2[3];
  emxInit_char_T(&iname, 2);
  if (fname->size[1] == 0) {
    i19 = iname->size[0] * iname->size[1];
    iname->size[0] = 1;
    iname->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)iname, i19, (int)sizeof(char));
  } else {
    in_idx = 0;
    copyfrom = 1;
    out_idx = 0;
    b_in_idx = 0;
    b_copyfrom = 1;
    out_len = 0;
    while (b_in_idx < fname->size[1]) {
      b_in_idx++;
      patt_idx = 1;
      tmp_in_idx = b_in_idx;
      while ((patt_idx <= 4) && (tmp_in_idx <= fname->size[1]) && (fname->
              data[tmp_in_idx - 1] == cv16[patt_idx - 1])) {
        tmp_in_idx++;
        patt_idx++;
      }

      if (patt_idx > 4) {
        out_len += 4;
        b_copyfrom = tmp_in_idx;
      } else {
        if (b_in_idx >= b_copyfrom) {
          out_len++;
        }
      }
    }

    i19 = iname->size[0] * iname->size[1];
    iname->size[0] = 1;
    iname->size[1] = out_len;
    emxEnsureCapacity((emxArray__common *)iname, i19, (int)sizeof(char));
    while (in_idx < fname->size[1]) {
      in_idx++;
      patt_idx = 1;
      tmp_in_idx = in_idx;
      while ((patt_idx <= 4) && (tmp_in_idx <= fname->size[1]) && (fname->
              data[tmp_in_idx - 1] == cv16[patt_idx - 1])) {
        tmp_in_idx++;
        patt_idx++;
      }

      if (patt_idx > 4) {
        for (b_in_idx = 0; b_in_idx < 4; b_in_idx++) {
          iname->data[out_idx] = cv17[b_in_idx];
          out_idx++;
        }

        copyfrom = tmp_in_idx;
      } else {
        if (in_idx >= copyfrom) {
          iname->data[out_idx] = fname->data[in_idx - 1];
          out_idx++;
        }
      }
    }
  }

  emxInit_int16_T(&t1_dim, 2);
  b_emxInit_int16_T(&r1, 1);
  b_emxInit_int32_T(&fp, 1);
  emxInit_uint8_T(&b_fp, 1);
  emxInit_uint8_T(&c_fp, 1);
  b_emxInit_int32_T(&d_fp, 1);
  b_emxInit_int16_T(&e_fp, 1);
  emxInit_uint8_T(&f_fp, 1);
  emxInit_uint8_T(&g_fp, 1);

  //  Open it if possible
  h_fp = f_fopen(fname);

  //  Read the fields
  e_fread(h_fp, fp);
  f_fread(h_fp, b_fp);
  g_fread(h_fp, c_fp);
  e_fread(h_fp, d_fp);
  h_fread(h_fp, e_fp);
  i_fread(h_fp, f_fp);
  i_fread(h_fp, g_fp);
  j_fread(h_fp, r1);
  i19 = t1_dim->size[0] * t1_dim->size[1];
  t1_dim->size[0] = 1;
  t1_dim->size[1] = r1->size[0];
  emxEnsureCapacity((emxArray__common *)t1_dim, i19, (int)sizeof(short));
  b_in_idx = r1->size[0];
  emxFree_uint8_T(&g_fp);
  emxFree_uint8_T(&f_fp);
  emxFree_int16_T(&e_fp);
  emxFree_int32_T(&d_fp);
  emxFree_uint8_T(&c_fp);
  emxFree_uint8_T(&b_fp);
  emxFree_int32_T(&fp);
  for (i19 = 0; i19 < b_in_idx; i19++) {
    t1_dim->data[t1_dim->size[0] * i19] = r1->data[i19];
  }

  emxFree_int16_T(&r1);
  emxInit_real32_T(&t1_pixdim, 2);
  b_emxInit_real32_T(&r2, 1);
  b_emxInit_real32_T(&i_fp, 1);
  b_emxInit_real32_T(&j_fp, 1);
  b_emxInit_real32_T(&k_fp, 1);
  b_emxInit_int16_T(&l_fp, 1);
  b_emxInit_int16_T(&m_fp, 1);
  b_emxInit_int16_T(&n_fp, 1);
  b_emxInit_int16_T(&o_fp, 1);
  k_fread(h_fp, i_fp);
  k_fread(h_fp, j_fp);
  k_fread(h_fp, k_fp);
  h_fread(h_fp, l_fp);
  h_fread(h_fp, m_fp);
  h_fread(h_fp, n_fp);
  h_fread(h_fp, o_fp);
  l_fread(h_fp, r2);
  i19 = t1_pixdim->size[0] * t1_pixdim->size[1];
  t1_pixdim->size[0] = 1;
  t1_pixdim->size[1] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)t1_pixdim, i19, (int)sizeof(float));
  b_in_idx = r2->size[0];
  emxFree_int16_T(&o_fp);
  emxFree_int16_T(&n_fp);
  emxFree_int16_T(&m_fp);
  emxFree_int16_T(&l_fp);
  emxFree_real32_T(&k_fp);
  emxFree_real32_T(&j_fp);
  emxFree_real32_T(&i_fp);
  for (i19 = 0; i19 < b_in_idx; i19++) {
    t1_pixdim->data[t1_pixdim->size[0] * i19] = r2->data[i19];
  }

  emxInit_real32_T(&t1_quatern_b, 2);
  b_emxInit_real32_T(&p_fp, 1);
  b_emxInit_real32_T(&q_fp, 1);
  b_emxInit_real32_T(&r_fp, 1);
  b_emxInit_int16_T(&s_fp, 1);
  emxInit_uint8_T(&t_fp, 1);
  emxInit_uint8_T(&u_fp, 1);
  b_emxInit_real32_T(&v_fp, 1);
  b_emxInit_real32_T(&w_fp, 1);
  b_emxInit_real32_T(&x_fp, 1);
  b_emxInit_real32_T(&y_fp, 1);
  b_emxInit_int32_T(&ab_fp, 1);
  b_emxInit_int32_T(&bb_fp, 1);
  emxInit_uint8_T(&cb_fp, 1);
  emxInit_uint8_T(&db_fp, 1);
  b_emxInit_int16_T(&eb_fp, 1);
  b_emxInit_int16_T(&fb_fp, 1);
  k_fread(h_fp, p_fp);
  k_fread(h_fp, q_fp);
  k_fread(h_fp, r_fp);
  h_fread(h_fp, s_fp);
  i_fread(h_fp, t_fp);
  i_fread(h_fp, u_fp);
  k_fread(h_fp, v_fp);
  k_fread(h_fp, w_fp);
  k_fread(h_fp, x_fp);
  k_fread(h_fp, y_fp);
  e_fread(h_fp, ab_fp);
  e_fread(h_fp, bb_fp);
  m_fread(h_fp, cb_fp);
  n_fread(h_fp, db_fp);
  h_fread(h_fp, eb_fp);
  h_fread(h_fp, fb_fp);
  k_fread(h_fp, r2);
  i19 = t1_quatern_b->size[0] * t1_quatern_b->size[1];
  t1_quatern_b->size[0] = 1;
  t1_quatern_b->size[1] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)t1_quatern_b, i19, (int)sizeof(float));
  b_in_idx = r2->size[0];
  emxFree_int16_T(&fb_fp);
  emxFree_int16_T(&eb_fp);
  emxFree_uint8_T(&db_fp);
  emxFree_uint8_T(&cb_fp);
  emxFree_int32_T(&bb_fp);
  emxFree_int32_T(&ab_fp);
  emxFree_real32_T(&y_fp);
  emxFree_real32_T(&x_fp);
  emxFree_real32_T(&w_fp);
  emxFree_real32_T(&v_fp);
  emxFree_uint8_T(&u_fp);
  emxFree_uint8_T(&t_fp);
  emxFree_int16_T(&s_fp);
  emxFree_real32_T(&r_fp);
  emxFree_real32_T(&q_fp);
  emxFree_real32_T(&p_fp);
  for (i19 = 0; i19 < b_in_idx; i19++) {
    t1_quatern_b->data[t1_quatern_b->size[0] * i19] = r2->data[i19];
  }

  emxInit_real32_T(&t1_quatern_c, 2);
  k_fread(h_fp, r2);
  i19 = t1_quatern_c->size[0] * t1_quatern_c->size[1];
  t1_quatern_c->size[0] = 1;
  t1_quatern_c->size[1] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)t1_quatern_c, i19, (int)sizeof(float));
  b_in_idx = r2->size[0];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    t1_quatern_c->data[t1_quatern_c->size[0] * i19] = r2->data[i19];
  }

  emxInit_real32_T(&t1_quatern_d, 2);
  k_fread(h_fp, r2);
  i19 = t1_quatern_d->size[0] * t1_quatern_d->size[1];
  t1_quatern_d->size[0] = 1;
  t1_quatern_d->size[1] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)t1_quatern_d, i19, (int)sizeof(float));
  b_in_idx = r2->size[0];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    t1_quatern_d->data[t1_quatern_d->size[0] * i19] = r2->data[i19];
  }

  emxInit_real32_T(&t1_qoffset_x, 2);
  k_fread(h_fp, r2);
  i19 = t1_qoffset_x->size[0] * t1_qoffset_x->size[1];
  t1_qoffset_x->size[0] = 1;
  t1_qoffset_x->size[1] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)t1_qoffset_x, i19, (int)sizeof(float));
  b_in_idx = r2->size[0];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    t1_qoffset_x->data[t1_qoffset_x->size[0] * i19] = r2->data[i19];
  }

  emxInit_real32_T(&t1_qoffset_y, 2);
  k_fread(h_fp, r2);
  i19 = t1_qoffset_y->size[0] * t1_qoffset_y->size[1];
  t1_qoffset_y->size[0] = 1;
  t1_qoffset_y->size[1] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)t1_qoffset_y, i19, (int)sizeof(float));
  b_in_idx = r2->size[0];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    t1_qoffset_y->data[t1_qoffset_y->size[0] * i19] = r2->data[i19];
  }

  emxInit_real32_T(&t1_qoffset_z, 2);
  k_fread(h_fp, r2);
  i19 = t1_qoffset_z->size[0] * t1_qoffset_z->size[1];
  t1_qoffset_z->size[0] = 1;
  t1_qoffset_z->size[1] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)t1_qoffset_z, i19, (int)sizeof(float));
  b_in_idx = r2->size[0];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    t1_qoffset_z->data[t1_qoffset_z->size[0] * i19] = r2->data[i19];
  }

  emxFree_real32_T(&r2);
  b_emxInit_real_T(&dim, 2);
  b_emxInit_real32_T(&gb_fp, 1);
  b_emxInit_real32_T(&hb_fp, 1);
  b_emxInit_real32_T(&ib_fp, 1);
  emxInit_uint8_T(&jb_fp, 1);
  emxInit_uint8_T(&kb_fp, 1);
  o_fread(h_fp, gb_fp);
  o_fread(h_fp, hb_fp);
  o_fread(h_fp, ib_fp);
  p_fread(h_fp, jb_fp);
  q_fread(h_fp, kb_fp);
  b_fclose(h_fp);
  i19 = dim->size[0] * dim->size[1];
  dim->size[0] = 1;
  dim->size[1] = t1_dim->size[1];
  emxEnsureCapacity((emxArray__common *)dim, i19, (int)sizeof(double));
  b_in_idx = t1_dim->size[0] * t1_dim->size[1];
  emxFree_uint8_T(&kb_fp);
  emxFree_uint8_T(&jb_fp);
  emxFree_real32_T(&ib_fp);
  emxFree_real32_T(&hb_fp);
  emxFree_real32_T(&gb_fp);
  for (i19 = 0; i19 < b_in_idx; i19++) {
    dim->data[i19] = t1_dim->data[i19];
  }

  emxFree_int16_T(&t1_dim);
  b_emxInit_real_T(&pixdim, 2);

  // Decode qform info from NIFTI-1 headers
  i19 = pixdim->size[0] * pixdim->size[1];
  pixdim->size[0] = 1;
  pixdim->size[1] = t1_pixdim->size[1];
  emxEnsureCapacity((emxArray__common *)pixdim, i19, (int)sizeof(double));
  b_in_idx = t1_pixdim->size[0] * t1_pixdim->size[1];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    pixdim->data[i19] = t1_pixdim->data[i19];
  }

  emxFree_real32_T(&t1_pixdim);
  b_emxInit_real_T(&b_t1_quatern_b, 2);

  //  Rotations from quaternions
  i19 = b_t1_quatern_b->size[0] * b_t1_quatern_b->size[1];
  b_t1_quatern_b->size[0] = 1;
  b_t1_quatern_b->size[1] = (t1_quatern_b->size[1] + t1_quatern_c->size[1]) +
    t1_quatern_d->size[1];
  emxEnsureCapacity((emxArray__common *)b_t1_quatern_b, i19, (int)sizeof(double));
  b_in_idx = t1_quatern_b->size[1];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    b_t1_quatern_b->data[b_t1_quatern_b->size[0] * i19] = t1_quatern_b->
      data[t1_quatern_b->size[0] * i19];
  }

  b_in_idx = t1_quatern_c->size[1];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    b_t1_quatern_b->data[b_t1_quatern_b->size[0] * (i19 + t1_quatern_b->size[1])]
      = t1_quatern_c->data[t1_quatern_c->size[0] * i19];
  }

  b_in_idx = t1_quatern_d->size[1];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    b_t1_quatern_b->data[b_t1_quatern_b->size[0] * ((i19 + t1_quatern_b->size[1])
      + t1_quatern_c->size[1])] = t1_quatern_d->data[t1_quatern_d->size[0] * i19];
  }

  emxFree_real32_T(&t1_quatern_d);
  emxFree_real32_T(&t1_quatern_c);
  emxFree_real32_T(&t1_quatern_b);
  emxInit_real_T(&img, 1);
  b_rot_matrix(b_t1_quatern_b, R);

  //  Translations
  i19 = img->size[0];
  img->size[0] = ((t1_qoffset_x->size[1] + t1_qoffset_y->size[1]) +
                  t1_qoffset_z->size[1]) + 1;
  emxEnsureCapacity((emxArray__common *)img, i19, (int)sizeof(double));
  b_in_idx = t1_qoffset_x->size[1];
  emxFree_real_T(&b_t1_quatern_b);
  for (i19 = 0; i19 < b_in_idx; i19++) {
    img->data[i19] = t1_qoffset_x->data[t1_qoffset_x->size[0] * i19];
  }

  b_in_idx = t1_qoffset_y->size[1];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    img->data[i19 + t1_qoffset_x->size[1]] = t1_qoffset_y->data
      [t1_qoffset_y->size[0] * i19];
  }

  b_in_idx = t1_qoffset_z->size[1];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    img->data[(i19 + t1_qoffset_x->size[1]) + t1_qoffset_y->size[1]] =
      t1_qoffset_z->data[t1_qoffset_z->size[0] * i19];
  }

  b_emxInit_real_T(&b_img, 2);
  img->data[(t1_qoffset_x->size[1] + t1_qoffset_y->size[1]) + t1_qoffset_z->
    size[1]] = 1.0;
  eye(dv1);
  b_in_idx = img->size[0];
  i19 = b_img->size[0] * b_img->size[1];
  b_img->size[0] = b_in_idx;
  b_img->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)b_img, i19, (int)sizeof(double));
  emxFree_real32_T(&t1_qoffset_z);
  emxFree_real32_T(&t1_qoffset_y);
  emxFree_real32_T(&t1_qoffset_x);
  for (i19 = 0; i19 < b_in_idx; i19++) {
    b_img->data[i19] = img->data[i19];
  }

  for (i19 = 0; i19 < 3; i19++) {
    for (copyfrom = 0; copyfrom < 4; copyfrom++) {
      T[copyfrom + (i19 << 2)] = dv1[copyfrom + (i19 << 2)];
    }
  }

  for (i19 = 0; i19 < 4; i19++) {
    T[12 + i19] = b_img->data[i19];
  }

  emxFree_real_T(&b_img);

  //  Zooms.  Note that flips are derived from the first
  //  element of pixdim, which is normally unused.
  if (dim->data[0] <= 3.0) {
    b_in_idx = (int)dim->data[0] + 1;
  } else {
    b_in_idx = 4;
  }

  if (2 > b_in_idx) {
    i19 = 1;
    copyfrom = 1;
  } else {
    i19 = 2;
    copyfrom = b_in_idx + 1;
  }

  b_emxInit_real_T(&Z, 2);
  b_copyfrom = 5 - b_in_idx;
  patt_idx = Z->size[0] * Z->size[1];
  Z->size[0] = 1;
  Z->size[1] = (copyfrom - i19) + b_copyfrom;
  emxEnsureCapacity((emxArray__common *)Z, patt_idx, (int)sizeof(double));
  b_in_idx = copyfrom - i19;
  for (patt_idx = 0; patt_idx < b_in_idx; patt_idx++) {
    Z->data[Z->size[0] * patt_idx] = pixdim->data[(i19 + patt_idx) - 1];
  }

  for (patt_idx = 0; patt_idx < b_copyfrom; patt_idx++) {
    Z->data[Z->size[0] * ((patt_idx + copyfrom) - i19)] = 1.0;
  }

  emxInit_boolean_T(&b_Z, 2);
  i19 = b_Z->size[0] * b_Z->size[1];
  b_Z->size[0] = 1;
  b_Z->size[1] = Z->size[1];
  emxEnsureCapacity((emxArray__common *)b_Z, i19, (int)sizeof(boolean_T));
  b_in_idx = Z->size[0] * Z->size[1];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    b_Z->data[i19] = (Z->data[i19] < 0.0);
  }

  emxInit_int32_T(&r3, 2);
  eml_li_find(b_Z, r3);
  b_in_idx = r3->size[0] * r3->size[1];
  emxFree_boolean_T(&b_Z);
  for (i19 = 0; i19 < b_in_idx; i19++) {
    Z->data[r3->data[i19] - 1] = 1.0;
  }

  emxFree_int32_T(&r3);
  if (pixdim->data[0] < 0.0) {
    Z->data[2] = -Z->data[2];
  }

  emxFree_real_T(&pixdim);
  b_emxInit_real_T(&c_Z, 2);
  diag(Z, c_Z);
  emxFree_real_T(&Z);
  for (i19 = 0; i19 < 4; i19++) {
    for (copyfrom = 0; copyfrom < 4; copyfrom++) {
      y[i19 + (copyfrom << 2)] = 0.0;
      for (patt_idx = 0; patt_idx < 4; patt_idx++) {
        y[i19 + (copyfrom << 2)] += T[i19 + (patt_idx << 2)] * R[patt_idx +
          (copyfrom << 2)];
      }
    }
  }

  b_emxInit_real_T(&M, 2);
  if (c_Z->size[0] == 1) {
    i19 = M->size[0] * M->size[1];
    M->size[0] = 4;
    M->size[1] = c_Z->size[1];
    emxEnsureCapacity((emxArray__common *)M, i19, (int)sizeof(double));
    for (i19 = 0; i19 < 4; i19++) {
      b_in_idx = c_Z->size[1];
      for (copyfrom = 0; copyfrom < b_in_idx; copyfrom++) {
        M->data[i19 + M->size[0] * copyfrom] = 0.0;
        for (patt_idx = 0; patt_idx < 4; patt_idx++) {
          M->data[i19 + M->size[0] * copyfrom] += y[i19 + (patt_idx << 2)] *
            c_Z->data[patt_idx + c_Z->size[0] * copyfrom];
        }
      }
    }
  } else {
    varargin_1_idx_1 = (unsigned int)c_Z->size[1];
    i19 = M->size[0] * M->size[1];
    M->size[0] = 4;
    emxEnsureCapacity((emxArray__common *)M, i19, (int)sizeof(double));
    i19 = M->size[0] * M->size[1];
    M->size[1] = (int)varargin_1_idx_1;
    emxEnsureCapacity((emxArray__common *)M, i19, (int)sizeof(double));
    b_in_idx = (int)varargin_1_idx_1 << 2;
    for (i19 = 0; i19 < b_in_idx; i19++) {
      M->data[i19] = 0.0;
    }

    if (c_Z->size[1] == 0) {
    } else {
      b_in_idx = (c_Z->size[1] - 1) << 2;
      for (copyfrom = 0; copyfrom <= b_in_idx; copyfrom += 4) {
        for (in_idx = copyfrom + 1; in_idx <= copyfrom + 4; in_idx++) {
          M->data[in_idx - 1] = 0.0;
        }
      }

      b_copyfrom = 0;
      for (copyfrom = 0; copyfrom <= b_in_idx; copyfrom += 4) {
        patt_idx = 0;
        for (out_len = b_copyfrom; out_len + 1 <= b_copyfrom + 4; out_len++) {
          if (c_Z->data[out_len] != 0.0) {
            tmp_in_idx = patt_idx;
            for (in_idx = copyfrom; in_idx + 1 <= copyfrom + 4; in_idx++) {
              tmp_in_idx++;
              M->data[in_idx] += c_Z->data[out_len] * y[tmp_in_idx - 1];
            }
          }

          patt_idx += 4;
        }

        b_copyfrom += 4;
      }
    }
  }

  emxFree_real_T(&c_Z);

  //  Convert from first voxel at [1,1,1]
  //  to first voxel at [0,0,0]
  if (M->size[1] == 1) {
    for (i19 = 0; i19 < 4; i19++) {
      for (copyfrom = 0; copyfrom < 4; copyfrom++) {
        R[i19 + (copyfrom << 2)] = 0.0;
        for (patt_idx = 0; patt_idx < 4; patt_idx++) {
          R[i19 + (copyfrom << 2)] += M->data[i19 + (patt_idx << 2)] * (double)
            b[patt_idx + (copyfrom << 2)];
        }
      }
    }
  } else {
    memset(&R[0], 0, sizeof(double) << 4);
    for (copyfrom = 0; copyfrom < 14; copyfrom += 4) {
      for (in_idx = copyfrom; in_idx + 1 <= copyfrom + 4; in_idx++) {
        R[in_idx] = 0.0;
      }
    }

    b_copyfrom = 0;
    for (copyfrom = 0; copyfrom < 14; copyfrom += 4) {
      patt_idx = 0;
      i19 = b_copyfrom + M->size[1];
      for (out_len = b_copyfrom; out_len + 1 <= i19; out_len++) {
        if (iv1[out_len] != 0) {
          tmp_in_idx = patt_idx;
          for (in_idx = copyfrom; in_idx + 1 <= copyfrom + 4; in_idx++) {
            tmp_in_idx++;
            R[in_idx] += (double)iv1[out_len] * M->data[tmp_in_idx - 1];
          }
        }

        patt_idx += 4;
      }

      b_copyfrom += M->size[1];
    }
  }

  emxFree_real_T(&M);

  // read volumes
  h_fp = g_fopen(iname);

  // fseek(fid, 0, 'bof');
  img_siz = dim->data[1];
  emxFree_char_T(&iname);
  for (b_in_idx = 0; b_in_idx < 6; b_in_idx++) {
    img_siz *= dim->data[b_in_idx + 2];
  }

  r_fread(h_fp, img_siz, img);

  // 'int16'
  for (i19 = 0; i19 < 7; i19++) {
    varargin_1[i19] = (short)dim->data[i19 + 1];
  }

  emxFree_real_T(&dim);
  for (i19 = 0; i19 < 7; i19++) {
    uv0[i19] = (unsigned short)varargin_1[i19];
  }

  d_emxInit_real_T(&c_img, 7);
  i19 = c_img->size[0] * c_img->size[1] * c_img->size[2] * c_img->size[3] *
    c_img->size[4] * c_img->size[5] * c_img->size[6];
  c_img->size[0] = uv0[0];
  c_img->size[1] = uv0[1];
  c_img->size[2] = uv0[2];
  c_img->size[3] = uv0[3];
  c_img->size[4] = uv0[4];
  c_img->size[5] = uv0[5];
  c_img->size[6] = uv0[6];
  emxEnsureCapacity((emxArray__common *)c_img, i19, (int)sizeof(double));
  for (b_in_idx = 0; b_in_idx + 1 <= img->size[0]; b_in_idx++) {
    c_img->data[b_in_idx] = img->data[b_in_idx];
  }

  emxFree_real_T(&img);

  //  set up the axes of the volume in voxel coordinates
  mri->dim[0] = (unsigned short)c_img->size[0];
  mri->dim[1] = (unsigned short)c_img->size[1];
  mri->dim[2] = (unsigned short)c_img->size[2];
  i19 = mri->anatomy->size[0] * mri->anatomy->size[1] * mri->anatomy->size[2] *
    mri->anatomy->size[3] * mri->anatomy->size[4] * mri->anatomy->size[5] *
    mri->anatomy->size[6];
  mri->anatomy->size[0] = c_img->size[0];
  mri->anatomy->size[1] = c_img->size[1];
  mri->anatomy->size[2] = c_img->size[2];
  mri->anatomy->size[3] = c_img->size[3];
  mri->anatomy->size[4] = c_img->size[4];
  mri->anatomy->size[5] = c_img->size[5];
  mri->anatomy->size[6] = c_img->size[6];
  emxEnsureCapacity((emxArray__common *)mri->anatomy, i19, (int)sizeof(double));
  b_in_idx = c_img->size[0] * c_img->size[1] * c_img->size[2] * c_img->size[3] *
    c_img->size[4] * c_img->size[5] * c_img->size[6];
  for (i19 = 0; i19 < b_in_idx; i19++) {
    mri->anatomy->data[i19] = c_img->data[i19];
  }

  emxFree_real_T(&c_img);
  for (i19 = 0; i19 < 16; i19++) {
    mri->transform[i19] = R[i19];
  }

  for (i19 = 0; i19 < 8; i19++) {
    mri->coordsys[i19] = cv18[i19];
  }

  //  try to determine the units of the coordinate system
  corner_points(mri->dim, R, pos_head);
  idrange(pos_head, dv2);
  estimate_units(norm(dv2), mri->unit);
}

