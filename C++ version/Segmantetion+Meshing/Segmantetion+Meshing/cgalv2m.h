
#ifndef __CGALV2M_H__
#define __CGALV2M_H__

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
extern void cgalv2m(const unsigned char vol[8530021], double opt, double maxvol,
                    double node[100564], double elem[565750], double face[362368]);

#endif

