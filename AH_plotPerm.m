function [analysisStruct,sigMaskInterp] = AH_plotPerm(data_SesbyPix,contrastLogic,varargin)
% Usage example see CSRTT_AnimalGroup_FC_CGC.m
% AH 5/26/2021

if ~isempty(varargin)
    permutationOptions = varargin{1};
else
    permutationOptions = struct([]);
end
% Freq downsample ratio
if isfield(permutationOptions,'fdsRatio')
    fdsRatio = permutationOptions.fdsRatio;
else
    fdsRatio = 1;
end
% Time downsample  
if isfield(permutationOptions,'tdsRatio')
    tdsRatio = permutationOptions.tdsRatio;
else
    tdsRatio = 1;
end
% tvec
if isfield(permutationOptions,'tvec')
    tvec = permutationOptions.tvec;
end
% baseTwin
if isfield(permutationOptions,'baseTwin')
    baseTwin = permutationOptions.baseTwin;
end

if isfield(permutationOptions,'sigOptions')
    sigOptions = permutationOptions.sigOptions;
else
    sigOptions = struct('onlyPos',0,'thresholdType','size'); % if 0, do both pos and neg
end
%% Compute permutation and mask
if numel(contrastLogic) == 1 % 1 cond
    tmp = data_SesbyPix(:,1:fdsRatio:end,1:tdsRatio:end);
    mat = shiftdim(tmp,1);% switch dimension freq x t x nSes (checked correct)
else
    tmp = data_SesbyPix(:,:,1:fdsRatio:end,1:tdsRatio:end);
    mat = shiftdim(tmp,2);% switch dimension freq x t x nSes x conds (checked correct)
end
[analysisStruct] = permutation2d_AH(mat,contrastLogic,permutationOptions);

if strcmp(sigOptions.thresholdType,'size')
    sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
        analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
        analysisStruct.permutation.sig.minSize,'apply','size',analysisStruct.permutation.sig.size,sigOptions);
elseif strcmp(sigOptions.thresholdType, 'mass')
    sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
        analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
        analysisStruct.permutation.sig.minSize,'apply','mass',analysisStruct.permutation.sig.mass,sigOptions);
end 

% Create query grid
if tdsRatio == 1 && fdsRatio == 1 % no downsample
    sigMaskInterp = sigMask; 
else
    if numel(size(tmp)) == 3
        [Xq,Yq] = ndgrid(1:1/tdsRatio:size(tmp,3),1/fdsRatio:1/fdsRatio:size(tmp,2));
        sigMaskInterp = interp2(1:size(tmp,3),1:size(tmp,2),sigMask,Xq,Yq,'linear')';  %X=1:n and Y=1:m, where [m,n] = size(V)
    elseif numel(size(tmp)) == 4
        [Xq,Yq] = ndgrid(1:1/tdsRatio:size(tmp,4),1/fdsRatio:1/fdsRatio:size(tmp,3)); 
        sigMaskInterp = interp2(1:size(tmp,4),1:size(tmp,3),sigMask,Xq,Yq,'linear')';  %X=1:n and Y=1:m, where [m,n] = size(V)
    end
end

% Based on permutation result, get a smallish mask
sigMaskInterp(sigMaskInterp<0.5)=0; 
sigMaskInterp(sigMaskInterp>=0.5)=1; 
if exist('baseTwin') && exist('tvec')
    try % For stim alignment, baseTwin is outside the range of matrix, then skip this
    % set everything before baseline to be 0
    sigMaskInterp(:,tvec<baseTwin(2)) = 0; 
    catch
    end
end
end