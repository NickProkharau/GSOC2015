function [node,elem,face]=cgalv2m(vol,opt,maxvol)
%
% convert a binary (or multi-valued) volume to tetrahedral mesh

% input:
%	 vol: a volumetric binary image
%	 ix,iy,iz: subvolume selection indices in x,y,z directions
%	 opt: parameters for CGAL mesher, if opt is a structure, then
%	     opt.radbound: defines the maximum surface element size
%	     opt.angbound: defines the miminum angle of a surface triangle
%	     opt.distbound: defines the maximum distance between the 
%		 center of the surface bounding circle and center of the 
%		 element bounding sphere
%	     opt.reratio:  maximum radius-edge ratio
%	     if opt is a scalar, it only specifies radbound.
%	 maxvol: target maximum tetrahedral elem volume
%
% output:
%	 node: output, node coordinates of the tetrahedral mesh
%	 elem: output, element list of the tetrahedral mesh, the last 
%	      column is the region id
%	 face: output, mesh surface element list of the tetrahedral mesh
%	      the last column denotes the boundary ID
%	      note: each triangle will appear twice in the face list with each
%		    one attaches to each side of the interface. one can remove
%		    the redundant triangles by unique(face(:,1:3),'rows')
%

fprintf(1,'creating surface and tetrahedral mesh from a multi-domain volume ...\n');
exesuff='.exe';


ang=30;
ssize=6;
approx=0.5;
reratio=3;

if(~isstruct(opt))
	ssize=opt;
end

saveinr(vol,mwpath('pre_cgalmesh.inr'));
fname=mwpath('post_cgalmesh.mesh');
if(exist(fname)) 
	delete(fname); 
end;

randseed=hex2dec('623F9A9E'); % "U+623F U+9A9E"

format_maxvol='%f';
mcpaht=[fileparts(which(mfilename)) '\' 'bin' '\' ('cgalmesh')];
cmd=sprintf(['"%s%s" "%s" "%s" %f %f %f %f ' format_maxvol ' %d'],mcpaht,exesuff,...
    mwpath('pre_cgalmesh.inr'),mwpath('post_cgalmesh.mesh'),ang,ssize,...
    approx,reratio,maxvol,randseed);
system(cmd);

[node,elem,face]=readmedit(mwpath('post_cgalmesh.mesh'));

fprintf(1,'node number:\t%d\ntriangles:\t%d\ntetrahedra:\t%d\nregions:\t%d\n',...
    size(node,1),size(face,1),size(elem,1),length(unique(elem(:,end))));
fprintf(1,'surface and volume meshes complete\n');

