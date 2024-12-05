function [cluster_members,visited] = check_neighbors(cluster_members, t_values, chan_hood, thresh, visited, sign)
    % Recursively check neighbors of the current cluster members
    % Initialize a stack to avoid repeated recursion (iterative approach)
    stack = cluster_members;
    
    while ~isempty(stack)
        % Pop the current channel from the stack
        ch = stack(end);
        stack(end) = [];
        
        % Find neighbors of the current channel
        neighbors = find(chan_hood(ch, :) == 1);
        
        for i = 1:length(neighbors)
            if ~visited(neighbors(i))
                if (sign == 1 && t_values(neighbors(i)) >= thresh) || (sign == -1 && t_values(neighbors(i)) <= thresh)
                    % Add the neighbor to the cluster
                    cluster_members = [cluster_members, neighbors(i)];
                    visited(neighbors(i)) = 1;
                    
                    % Add the neighbor to the stack to check its neighbors
                    stack = [stack, neighbors(i)];
                end
            end
        end
    end
end