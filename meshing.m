% sample script to create volumetric mesh from 
% multiple levelsets of a binary segmented head image.
restoredefaultpath;
addpath ('mesh');
% load segmented brain images 
load segmentedmri.mat

head = segmentedmri.scalp;
brain=segmentedmri.brain;
scull=segmentedmri.skull;

% fill holes in the head image and create the canonical binary volume
fprintf(1,'filling holes in the volumetric images...\n');
tic
maxgap=10;
cleanimg = imclose(logical(head>0),strel(ones(maxgap,maxgap,maxgap)));
cleanimg=imfill(cleanimg,'holes');

cleanbrain = imclose(logical(brain>0),strel(ones(maxgap,maxgap,maxgap)));
cleanbrain=imfill(cleanbrain,'holes');

cleanscull = imclose(logical(scull>0),strel(ones(maxgap,maxgap,maxgap)));
cleanscull=imfill(cleanscull,'holes');

toc

% add brain image as additional segment
cleanimgfull=cleanimg+(cleanscull>0)+(cleanbrain>0);

clear cleanbrain;
clear cleanscull;
clear head;
clear brain;
clear scull;
clear segmentedmri;

% create volumetric tetrahedral mesh from the two-layer 3D images
clear opt;

opt(1).radbound=4; % head surface element size bound
opt(2).radbound=2; % brain surface element size bound
opt(3).radbound=2; % skull surface element size bound

%convert image to uint8
img=uint8(cleanimgfull);

%meshing
tic
vol=img(1:size(cleanimg,1),1:size(cleanimg,2),1:size(cleanimg,3));
[node, elem,face]=cgalv2m(img,opt,100);		
toc

h=slice((cleanimgfull),[],[120],[120]);
set(h,'linestyle','none')
figure;

plotmesh(node(:,[2 1 3]),face,'facealpha',0.7,'linestyle','none');
z0=mean(node(:,3));

plane=[min(node(:,1)) min(node(:,2)) z0
       min(node(:,1)) max(node(:,2)) z0
       max(node(:,1)) min(node(:,2)) z0];

% run qmeshcut to get the cross-section information at z=mean(node(:,1))
% use the x-coordinates as the nodal values
figure;
hsurf=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),'facecolor','none');

[cutpos,cutvalue,facedata]=qmeshcut(elem(:,1:4),node(:,1:3),node(:,1),plane);

figure;
hsurf=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),'facecolor','none');
hold on;
 patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','interp');


