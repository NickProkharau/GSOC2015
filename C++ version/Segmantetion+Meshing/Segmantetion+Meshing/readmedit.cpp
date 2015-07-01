
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
#include "fclose.h"
#include "Segmentation1_emxutil.h"
#include "fread.h"
#include "fopen.h"
#include <stdio.h>


// Function Definitions

//
// [node,elem,face]=readmedit(filename)
//
//  read Medit mesh format
// Arguments    : double node[100564]
//                double elem[565750]
//                double face[362368]
// Return Type  : void
//
void b_readmedit(double node[100564], double elem[565750], double face[362368])
{
  emxArray_real_T *b_node;
  double fid;
  int k;
  emxInit_real_T(&b_node, 1);

  
  fid = c_fopen();
  b_fread(fid, b_node);

  // 'int16'
  for (k = 0; k + 1 <= b_node->size[0]; k++) {
    node[k] = b_node->data[k];
  }

  b_fclose(fid);
  fid = d_fopen();
  c_fread(fid, b_node);

  // 'int16'
  for (k = 0; k + 1 <= b_node->size[0]; k++) {
    elem[k] = b_node->data[k];
  }

  b_fclose(fid);
  fid = e_fopen();
  d_fread(fid, b_node);

  // 'int16'
  for (k = 0; k + 1 <= b_node->size[0]; k++) {
    face[k] = b_node->data[k];
  }

  emxFree_real_T(&b_node);
  b_fclose(fid);
  for (k = 0; k < 100564; k++) {
    node[k] *= 112.0;
  }

  //  fid=fopen(filename,'r');
  //  while(1)
  //      key=fscanf(fid,'%s',1);
  //      if(strcmp(key,'End')) break; end
  //      val=fscanf(fid,'%d',1);
  //      if(strcmp(key,'Vertices'))
  //          node=fscanf(fid,'%f',4*val);
  //          node=reshape(node,[4 val])';
  //      elseif(strcmp(key,'Triangles'))
  //          face=fscanf(fid,'%d',4*val);
  //          face=reshape(face,[4 val])';
  //      elseif(strcmp(key,'Tetrahedra'))
  //          elem=fscanf(fid,'%d',5*val);
  //          elem=reshape(elem,[5 val])';
  //      end
  //  end
  //
  //  fclose(fid);
}

//
// [node,elem,face]=readmedit(filename)
//
//  read Medit mesh format
// Arguments    : const unsigned char vol[8530021]
//                double node[100564]
//                double elem[565750]
//                double face[362368]
// Return Type  : void
//
void readmedit(const unsigned char vol[8530021], double node[100564], double
               elem[565750], double face[362368])
{
  emxArray_real_T *b_node;
  double fid;
  int k;
  emxInit_real_T(&b_node, 1);

  //
  //  input:
  //     fname: name of the medit data file
  //
  //  output:
  //     node: node coordinates of the mesh
  //     elem: list of elements of the mesh	
  //     face: list of surface triangles of the mesh	
  //
  //  node=[];
  //  elem=[];
  //  face=[];
  //
  fid = c_fopen();
  b_fread(fid, b_node);

  // 'int16'
  for (k = 0; k + 1 <= b_node->size[0]; k++) {
    node[k] = b_node->data[k];
  }

  b_fclose(fid);
  fid = d_fopen();
  c_fread(fid, b_node);

  // 'int16'
  for (k = 0; k + 1 <= b_node->size[0]; k++) {
    elem[k] = b_node->data[k];
  }

  b_fclose(fid);
  fid = e_fopen();
  d_fread(fid, b_node);

  // 'int16'
  for (k = 0; k + 1 <= b_node->size[0]; k++) {
    face[k] = b_node->data[k];
  }

  emxFree_real_T(&b_node);
  b_fclose(fid);
  for (k = 0; k < 100564; k++) {
    node[k] *= (double)vol[0];
  }
  
}
