close all;clear all;clc;
%%
npoints = 5000;
nquery = 500;
ndims = 2;

imgsize = 640;

points_truth       = rand(npoints, ndims)*imgsize;
points_query_truth = points_truth(randsample(npoints,nquery),:) + 0*randn(nquery, ndims);

%% truth
binSize   = npoints;
distThres = 100.0;
tree      = BKDTree_build_mex(points_truth,binSize,distThres);

tic;
ref_truth = BKDTree_nearest_neighbor_mex( tree, points_query_truth);
BKDTree_delete_mex(tree);
toc;

%% bin size
ref_save = zeros(nquery, npoints);
diff_num_matching = zeros(npoints,1);
num_monte_carlo = 10;
distThres = 100.0;
dt_save = zeros(30,npoints);
diff_save = zeros(nquery,npoints,num_monte_carlo);
for j = 1:num_monte_carlo
    points       = rand(npoints, ndims)*imgsize;
    points_query = points(randsample(npoints, nquery),:) + 1.5*randn(nquery, ndims);
    tree         = BKDTree_build_mex(points, npoints, distThres);
    ref_truth    = BKDTree_nearest_neighbor_mex( tree, points_query);
    BKDTree_delete_mex(tree);
    
    for i = 1:40:npoints
        fprintf('# of monte carlo : %d, # of binsize : %d\n',j,i);
        binSize = i;
        tree      = BKDTree_build_mex(points, binSize, distThres);
        
        if(j == num_monte_carlo)
            for k=1:30
                tic;
                ref_temp  = BKDTree_nearest_neighbor_mex( tree, points_query );
                dt = toc;
                dt_save(k,i) = dt;
            end
            
        else
            ref_temp  = BKDTree_nearest_neighbor_mex( tree, points_query);
        end
        
        diff_save(:,i,j) = (ref_temp - ref_truth) ~= 0;
        
        BKDTree_delete_mex(tree);
    end
end
%%
for i = 1: npoints
    dt_mean(i) = mean(dt_save(:,i));
    dt_std(i)  = std(dt_save(:,i));
    diff_temp = diff_save(:,i,:);
    diff_num_matching(i) = sum(diff_temp(:)) / num_monte_carlo;
end

success_rate = 100 - diff_num_matching/nquery*100;
efficiency = success_rate ./ dt_mean.';

figure();
plot(100 - diff_num_matching/nquery*100,'k','linewidth',2);xlabel('bin size'); ylabel('matching rate [percent]');ylim([0,100]);
title('matching rate vs. bin size');

figure();
plot(dt_mean*1000);xlabel('bin size'); ylabel('time consumption per iteration [ms]');
title('time consumption vs. bin size');

figure();
efficiency(find(efficiency>1e10))=0;
plot(efficiency);xlabel('bin size');ylabel('efficiency');
title('efficiency vs. bin size');
%% draw points & lines
% figure();
% plot(points_truth(:,1),points_truth(:,2),'r*'); hold on; grid on;
% plot(points_query_truth(:,1),points_query_truth(:,2),'b.');
% for i = 1:length(ref_truth)
%     line([points_truth(ref_truth(i),1),points_query_truth(i,1)],[points_truth(ref_truth(i),2),points_query_truth(i,2)]);
% end
% % a=text(-20,-70,'bin size : %d\n','fontsize',12);
% 
% offsetsize = 100;
% xlim([-offsetsize,imgsize+offsetsize]);ylim([-offsetsize,imgsize+offsetsize]); legend('references','queries');
% 
