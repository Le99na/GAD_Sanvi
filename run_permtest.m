
dat_eng = load('halpha_dis_pwr_clst_frmt.mat');
dat_eng_3d = dat_eng.data;

dat_dis = load('halpha_eng_pwr_clst_frmt.mat');
dat_dis_3d = dat_dis.data(:,1:5851,:);

max_dist = 0.4

loadedData = load('chanlocs.mat');

chan_hood = spatial_neighbors(loadedData.chan_hood,max_dist);

thresh_p = 0.05
verblevel = 2
%seed_state = 3
n_perm = 5000
fwer = 0.05
tail = -1
freq_domain = 0

[pval, t_orig, clust_info, ~, est_alpha] = clust_perm2_dependent(eng,dis,chan_hood,n_perm,fwer,tail,thresh_p,verblevel,[],freq_domain)

% Save the variables to a .mat file
save('/work/users/m/a/magcam/channels_2nd/cluster_output/PLI_theta_eng_diseng.mat', 'pval', 't_orig', 'clust_info', 'est_alpha');