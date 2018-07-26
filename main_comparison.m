close all;clear all;clc;
%% parameter settings
npoints  = 500;
nquery   = 500;
ndims     = 5;

imgsize = 320;

binSize       = npoints;
dist_thres = 100.0;


points_truth       = rand(npoints, ndims)*imgsize;
points_query_truth = points_truth(randsample(npoints,nquery),:) + 1.5*randn(nquery, ndims);

%% truth

bkdtree_root = BKDTree_build_mex(points_truth,binSize,dist_thres);
ref_truth_bkdtree = BKDTree_nearest_neighbor_mex( bkdtree_root, points_query_truth);
BKDTree_delete_mex(bkdtree_root);

%% 1. BKDTree & KDTree : according to bin size
ref_save = zeros(nquery, npoints);
num_monte_carlo  = 30;
dist_thres = 100000000.0;

bin_gap  = 5;
bin_max = 200;
bin_index_max = bin_max/bin_gap;
bin_index = bin_gap:bin_gap:bin_max;

dt_bkd_save = zeros(30, bin_index_max);
dt_kd_save = zeros(30, bin_index_max);

diff_bkd_save = zeros(nquery, bin_index_max, num_monte_carlo);
diff_kd_save = zeros(nquery, bin_index_max, num_monte_carlo);

diff_bkd_num_matching = zeros(bin_index_max, 1);
diff_kd_num_matching = zeros(bin_index_max, 1);

for j = 1:num_monte_carlo
    points = rand(npoints, ndims)*imgsize;
    points_query = points(randsample(npoints, nquery),:) + 2.5*randn(nquery, ndims);
    
    bkdtree_root = BKDTree_build_mex(points, npoints, dist_thres);
    kdtree_root   = kdtree_build(points);
    
    ref_truth_bkdtree    = BKDTree_nearest_neighbor_mex( bkdtree_root, points_query);
    ref_truth_kdtree      = kdtree_nearest_neighbor( kdtree_root, points_query, 1e99);

    BKDTree_delete_mex(bkdtree_root);
    kdtree_delete(kdtree_root);
    
    for i = 1:bin_max/bin_gap
        fprintf('# of monte carlo : %d, # of binsize : %d\n',j,i);
        binSize = bin_index(i);
        kdtree_root   = kdtree_build(points);
        bkdtree_root = BKDTree_build_mex(points, binSize, dist_thres);
        
        if(j == num_monte_carlo)
            for k=1:100
                tic;
                ref_bkd_temp  = BKDTree_nearest_neighbor_mex( bkdtree_root, points_query );
                dt = toc;
                dt_bkd_save(k,i) = dt;
                
                tic;
                ref_kd_temp    = kdtree_nearest_neighbor( kdtree_root, points_query, dist_thres);
                dt = toc;
                dt_kd_save(k,i) = dt;
            end
        else
            ref_bkd_temp  = BKDTree_nearest_neighbor_mex( bkdtree_root, points_query);
            ref_kd_temp    = kdtree_nearest_neighbor( kdtree_root, points_query, dist_thres);
        end
        
        diff_bkd_save(:,i,j) = (ref_bkd_temp - ref_truth_bkdtree) ~= 0;
        diff_kd_save(:,i,j)   = (ref_kd_temp - ref_truth_kdtree) ~= 0;

        BKDTree_delete_mex(bkdtree_root);
        kdtree_delete(kdtree_root);
    end
end
%% graphs
for i = 1: bin_index_max
    dt_bkd_mean(i) = mean(dt_bkd_save(:,i));
    dt_bkd_std(i)     = std(dt_bkd_save(:,i));
    
    dt_kd_mean(i) = mean(dt_kd_save(:,i));
    dt_kd_std(i)     = std(dt_kd_save(:,i));
    
    diff_temp = diff_bkd_save(:,i,:);
    diff_bkd_num_matching(i) = sum(diff_temp(:)) / num_monte_carlo;
    
    diff_temp = diff_kd_save(:,i,:);
    diff_kd_num_matching(i) = sum(diff_temp(:)) / num_monte_carlo;
end

success_bkd_rate = 100 - diff_bkd_num_matching/nquery*100;
success_kd_rate   = 100 - diff_kd_num_matching/nquery*100;

figure();
plot(bin_index, success_bkd_rate,'k','linewidth',1);xlabel('bin size'); ylabel('matching rate [percent]');ylim([0,100]);
title('bkdtree - matching rate vs. bin size');

figure();
plot(bin_index,dt_bkd_mean*1000);xlabel('bin size'); ylabel('time consumption per iteration [ms]');
title('bkdtree - time consumption vs. bin size');

figure();
plot(bin_index, success_kd_rate,'k','linewidth',1);xlabel('bin size'); ylabel('matching rate [percent]');ylim([0,100]);
title('kdtree - matching rate vs. bin size');

figure();
plot(bin_index,dt_kd_mean*1000);xlabel('bin size'); ylabel('time consumption per iteration [ms]');
title('kdtree - time consumption vs. bin size');

figure();
plot(bin_index, dt_bkd_mean./dt_kd_mean,'k','linewidth',1); hold on; line([0,bin_index(end)],[1,1],'linestyle','--','linewidth',2);
xlim([0,bin_index(end)]); ylim([0,max(dt_bkd_mean./dt_kd_mean)]);
xlabel('bin size'); ylabel('time ratio (t_b_k_d / t_k_d)');
title('time consumption comparison');

figure();
plot(bin_index, success_bkd_rate.'./( dt_bkd_mean./dt_kd_mean));
%% draw points & lines
profile on;
points = rand(npoints, ndims)*imgsize;
points_query = points(randsample(npoints, nquery),:) + 2.5*randn(nquery, ndims);
kdtree_root = kdtree_build(points);
bkdtree_root = BKDTree_build_mex(points, 30,  999999999.999);
tic;
for i = 1:100
    ref_truth_bkdtree = kdtree_nearest_neighbor(kdtree_root, points_query, 999999999.9);
end
toc;

tic;
for i = 1:100
    ref_truth_bkdtree    = BKDTree_nearest_neighbor_mex( bkdtree_root, points_query );
end
toc;

kdtree_delete(kdtree_root);
figure();
plot(points(:,1),points(:,2),'r*'); hold on; grid on;
plot(points_query(:,1),points_query(:,2),'b.');
for i = 1:length(ref_truth_bkdtree)
    line([points(ref_truth_bkdtree(i),1),points_query(i,1)],[points(ref_truth_bkdtree(i),2),points_query(i,2)]);
end
% a=text(-20,-70,'bin size : %d\n','fontsize',12);

offsetsize = 100;
xlim([-offsetsize,imgsize+offsetsize]);ylim([-offsetsize,imgsize+offsetsize]); legend('references','queries');

profile viewer;