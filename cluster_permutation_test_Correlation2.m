function [r_values, p_values, sig_clusters] = cluster_permutation_test_Correlation2(data, ratings, n_permutations, alpha)
    % Input:
    %   data - PAC data (channels x subjects)
    %   ratings - questionnaire ratings (1 x subjects)
    %   n_permutations - number of permutations for cluster correction
    %   alpha - significance threshold (e.g., 0.05)
    %
    % Output:
    %   t_values - t-values for each channel
    %   p_values - p-values for each channel
    %   sig_clusters - significant clusters with corrected p-values

    % Load the neighborhood information (spatial location of channels)
    load('chanlocs_90.mat', 'data_chanlocs'); 
    chan_hood = spatial_neighbors2(data_chanlocs, 3.5);%0.4

    % Get the number of channels and subjects
    [n_channels, n_subjects] = size(data);

    % Check if the dimensions match
    if length(ratings) ~= n_subjects
        error('The dimensions of the data and ratings must match!');
    end

    % Calculate correlation and t-values for each channel between data and ratings
    t_values = zeros(n_channels, 1);
    p_values = zeros(n_channels, 1);
    r_values = zeros(n_channels, 1);
    for ch = 1:n_channels
        [r, p] = corr(data(ch, :)', ratings', 'Type', 'Spearman');
        t = sqrt((r^2 * (n_subjects - 2)) / (1 - r^2));  % Convert r to t
        t_values(ch) = t;
        p_values(ch) = p;
        r_values(ch) = r;
    end

    % Threshold the t-values based on the alpha level for a two-tailed test
    t_critical = tinv(1 - alpha / 2, n_subjects - 2);  % Critical t-value for two-tailed test

    % Find clusters of significant channels (positive and negative t-values)
    pos_clusters = find_clusters_mass(t_values, t_critical, chan_hood, 1);  % For positive clusters
    neg_clusters = find_clusters_mass(t_values, -t_critical, chan_hood, -1);  % For negative clusters

    % Combine positive and negative clusters
    all_clusters = [pos_clusters, neg_clusters];

    % Cluster-based permutation testing
    max_cluster_stats = zeros(n_permutations, 1);
    for perm = 1:n_permutations
        % Permute the ratings (shuffle across participants)
        perm_ratings = ratings(randperm(n_subjects));

        % Compute correlation values for permuted data
        perm_t_values = zeros(n_channels, 1);
        for ch = 1:n_channels
            r = corr(data(ch, :)', perm_ratings', 'Type', 'Spearman');
            perm_t_values(ch) = sqrt((r^2 * (n_subjects - 2)) / (1 - r^2));  % Convert r to t for permuted data
        end

        % Identify clusters in permuted data
        perm_pos_clusters = find_clusters_mass(perm_t_values, t_critical, chan_hood, 1);
        perm_neg_clusters = find_clusters_mass(perm_t_values, -t_critical, chan_hood, -1);

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