function bkdtree_compile()
% clc
localpath = fileparts(which('BKDTree_compile'));
fprintf(1,'Compiling BKDTree library [%s]...\n', localpath);

err = 0;
err = err | mex('-outdir',localpath, [localpath,'/BKDTree_build_mex.cpp']);
err = err | mex('-outdir',localpath, [localpath,'/BKDTree_nearest_neighbor_mex.cpp']);

if err ~= 0, 
   error('compile failed!'); 
else
   fprintf(1,'\bDone!\n');
end 