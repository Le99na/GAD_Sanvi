function clusters = find_clusters_mass(t_values, thresh, chan_hood, sign)
    % Function to find clusters of significant channels
    % sign: 1 for positive clusters, -1 for negative clusters
    
    n_channels = length(t_values);
    clusters = [];
    cluster_count = 0;
    
    % Identify channels above threshold
    if sign == 1
        sig_channels = find(t_values >= thresh);
    else
        sig_channels = find(t_values <= thresh);
    end
    
    % Keep track of visited channels
    visited = zeros(n_channels, 1);
    
    % Iterate over all significant channels
    for i = 1:length(sig_channels)
        ch = sig_channels(i);
        if ~visited(ch)
            % Start a new cluster
            cluster_count = cluster_count + 1;
            cluster_members = ch;
            visited(ch) = 1;
            
            % Check neighbors recursively
            [cluster_members, visited] = check_neighbors(cluster_members, t_values, chan_hood, thresh, visited, sign);
            
            % Save the cluster
            clusters(cluster_count).channels = cluster_members;
            clusters(cluster_count).cluster_mass = sum(abs(t_values(cluster_members))); % Changed to cluster mass
        end
    end
end