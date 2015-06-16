function tempname=mwpath(fname)

% get full temp-file name by prepend working-directory and current session name

username=getenv('UserName'); % for windows

username=['mesh-' username];

tdir=tempdir;
tdir=[tdir username filesep];
tempname=[tdir fname];
