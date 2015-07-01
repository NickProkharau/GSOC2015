Segmentation and meshing scripts by Prohorov N. A.

These scripts are written in MatLab for Windows.
MRI images of head mni_icbm152_t1_tal_nlin_asym_09c are used as Input files. 
mni_icbm152_t1_tal_nlin_asym_09c.hdr - input header file of MRI images.
mni_icbm152_t1_tal_nlin_asym_09c.img - file with MRI images.

There are 2 main scripts: segmentation.m and meshing.m

segmentation.m script reads input images in special form for further processing. Then It makes a segmentation of head volume into scalp, skull and brain and saves obtained segmentation in segmentationmri.mat file.

meshing.m script load segmentationmri.mat file and process it. First this script fills holes in each part of segmentation (scalp, skull and brain), combines scalp, skull and brain into one volume and meshes this volume with tetrahedral elements. Nodes, elements and faces are output of meshing.  As the result script makes diffrent visualizations of mesh using nodes, elements and faces abtained before.

A C++ version of these scripts had been added