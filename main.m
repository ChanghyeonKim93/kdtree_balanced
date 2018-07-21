close all;clear all;clc;
%%
npoints = 500;
nquery = 500;
ndims = 2;

imgsize = 640;
binSize = 10;
distThres = 100.0;

points = rand(npoints, ndims)*imgsize;
points_query = points(randsample(npoints,nquery),:) + 1.5*randn(nquery, ndims);

tree = BKDTree_build_mex(points,binSize,distThres);

tic;
for k = 1:1
ref_ind = BKDTree_nearest_neighbor_mex( tree, points_query);
end
BKDTree_delete_mex(tree);
toc;

figure();
plot(points(:,1),points(:,2),'r*'); hold on; grid on;
plot(points_query(:,1),points_query(:,2),'b.');
for i = 1:length(ref_ind)
    line([points(ref_ind(i),1),points_query(i,1)],[points(ref_ind(i),2),points_query(i,2)]);
end
a=text(-20,-70,'bin size : %d\n','fontsize',12);

offsetsize = 100;
xlim([-offsetsize,imgsize+offsetsize]);ylim([-offsetsize,imgsize+offsetsize]); legend('references','queries');

memory