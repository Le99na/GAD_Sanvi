function [t_values, p_values, sig_clusters] = cluster_permutation_test_oneSample(data, n_permutations, alpha)
    % Input:
    %   data - dataset (channels x subjects)
    %   n_permutations - number of permutations for cluster correction
    %   alpha - significance threshold (e.g., 0.05)
    %
    % Output:
    %   t_values - t-values from one-sample t-test (channels x 1)
    %   p_values - p-values from one-sample t-test (channels x 1)
    %   sig_clusters - significant clusters with corrected p-values

    % Load the neighborhood information (spatial location of channels)
    load('chanlocs.mat', 'chan_hood');  % Assuming chan_hood is (channels x channels)
    chan_hood = spatial_neighbors2(chan_hood, 0.4);%3.5
    
    % Get the number of channels and subjects
    [n_channels, n_subjects] = size(data);
    
    % One-sample t-test against zero for each channel
    [~, p_values, ~, stats] = ttest(data', 0, 'alpha', alpha);
    t_values = stats.tstat;  % Extract the t-values

    % Threshold the t-values based on the alpha value for a two-tailed test
    thresh = tinv(1 - alpha / 2, n_subjects - 1);
    
    % Find clusters of significant channels (positive and negative t-values)
    pos_clusters = find_clusters_mass(t_values, thresh, chan_hood, 1);  % For positive clusters
    neg_clusters = find_clusters_mass(t_values, -thresh, chan_hood, -1);  % For negative clusters
    
    % Combine positive and negative clusters
    all_clusters = [pos_clusters, neg_clusters];
    
    % Cluster-based permutation testing
    max_cluster_stats = zeros(n_permutations, 1);
    for perm = 1:n_permutations
        
        % Shuffle the signs within each channel to create a random distribution
        perm_data = data .* (2 * (rand(n_channels, n_subjects) > 0.5) - 1);
        
        % Compute t-values for permuted data
        [~, ~, ~, perm_stats] = ttest(perm_data', 0);
        perm_t_values = perm_stats.tstat;

        % Identify clusters in permuted data
        perm_pos_clusters = find_clusters_mass(perm_t_values, thresh, chan_hood, 1);
        perm_neg_clusters = find_clusters_mass(perm_t_values, -thresh, chan_hood, -1);

        % Get the largest cluster statistic from permuted data
        if ~isempty(perm_pos_clusters)
            max_cluster_stats(perm) = max([perm_pos_clusters.cluster_mass]);
        end
        if ~isempty(perm_neg_clusters)
            max_cluster_stats(perm) = max([max_cluster_stats(perm); max([perm_neg_clusters.cluster_mass])]);
        end
    end
    
    % Compute corrected p-values for the original clusters
    sig_clusters = struct('cluster_idx', [], 'p_val', [], 'channels', []);
    cluster_count = 0;  
    for clust_idx = 1:length(all_clusters)
        clust_stat = all_clusters(clust_idx).cluster_mass;
        p_val = mean(max_cluster_stats >= clust_stat);
        
        % Check if the cluster is significant
        if p_val < alpha
            cluster_count = cluster_count + 1;
            sig_clusters(cluster_count).cluster_idx = clust_idx;
            sig_clusters(cluster_count).p_val = p_val;
            
            % Store the channels that belong to this significant cluster
            sig_clusters(cluster_count).channels = all_clusters(clust_idx).channels;
        end
    end
end