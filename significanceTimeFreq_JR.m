function [varargout] = significanceTimeFreq_JR(tValues,pValues,alphaThreshold,minClusterSize,useFlag,varargin)
% Two use cases: to test for significance and to apply that threshold
if strcmp(useFlag,'test')
    % output is maxSize, maxT, maxMass, maxSize
    maxSize = 0;
    maxT = 0;
    maxMass = 0;
    inputCounter = 0;
elseif strcmp(useFlag,'apply')
    % output is a mask
    sigType = varargin{1};
    assert(~isempty(find(strcmpi(sigType,{'size','mass','k'}),1)));
    if strcmp(sigType,'k')
        inputCounter = 1;
    else
        sigValue = varargin{2};
        inputCounter = 2;
    end
    sigMask = zeros(size(tValues));
else
    error('Invalid use case');
end
% Optional input
if length(varargin) > inputCounter
    sigOptions = varargin{inputCounter+1};
else
    sigOptions = struct([]);
end
% Only positive
if isfield(sigOptions,'onlyPos')
    onlyPos = sigOptions.onlyPos;
else
    onlyPos = 0;
end

% The p-map is misleading because values can be positive and negatively
% significant but next to each other
% thresold data by alpha
pMap = double(pValues < alphaThreshold);
% Run for positive and negative clusters separately
pos_pMap = pMap .* double(tValues > 0);
neg_pMap = pMap .* double(tValues < 0);
% find uncorrected signifiance clusters
pos_pClustersCC = bwconncomp(pos_pMap); % grab connected pixels as a cluster
pos_pClusters = pos_pClustersCC.PixelIdxList;
neg_pClustersCC = bwconncomp(neg_pMap);
neg_pClusters = neg_pClustersCC.PixelIdxList;
if onlyPos
    % Only look at positive clusters - one-tailed
    pClusters = pos_pClusters;
else
    % Combine negative and positive clusters
    pClusters = [pos_pClusters neg_pClusters];
end
% loop through clusters
for clusterIdx = 1:length(pClusters)
    cluster = pClusters{clusterIdx};
    
    % t-values for this cluster
    tCluster = abs(tValues(cluster));
    
    % Calculate size, max, and mass
    clusterSize = length(cluster); % size of the cluster
    clusterMaxT = max(tCluster); % max t-value of the cluster
    clusterMass = mean(tCluster); % average t-value within the cluster
    
    % enforce minimum cluster size to avoid warp from spurious clusters
    if clusterSize > minClusterSize
        % Two use cases:
        if strcmp(useFlag,'test')
            % output is maxT, maxMass, maxSize
            if clusterSize > maxSize
                maxSize = clusterSize;
            end
            if clusterMaxT > maxT
                maxT = clusterMaxT;
            end
            if clusterMass > maxMass
                maxMass = clusterMass;
            end
        elseif strcmp(useFlag,'apply')
            switch sigType
                case 'mass'
                    % Check for signifiance
                    if clusterMass > sigValue
                        sigMask(cluster) = 1;
                    end
                case 'size'
                    if clusterSize > sigValue
                        sigMask(cluster) = 1;
                    end
                case 'k'
                    sigMask(cluster) = 1;
            end
        end
    end % min cluster size
end % loop clusters

% Ouput
if strcmp(useFlag,'test')
    % output is maxSize, maxT, and maxMass
    varargout{1} = maxSize;
    varargout{2} = maxT;
    varargout{3} = maxMass;
elseif strcmp(useFlag,'apply')
    % output is a mask
    varargout{1} = sigMask;
end
end % end of function