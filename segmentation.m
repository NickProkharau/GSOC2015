restoredefaultpath;
addpath ('spm8');
addpath ('segm');

% reading MRI images
mri = ft_read_mri('mni_icbm152_t1_tal_nlin_asym_09c.hdr');
mri.coordsys = 'neuromag';
disp(mri);  

% segmentation
cfg           = [];
cfg.output    = {'brain','skull','scalp'};
save mri;
segmentedmri  = ft_volumesegment(cfg, mri);

save segmentedmri segmentedmri

