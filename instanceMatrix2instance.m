% instanceMatrix2instance.m
%
%      usage: [instances stimVals] = instanceMatrix2instance(instanceMatrix,stimVals)
%         by: justin gardner
%       date: 09/12/16
%    purpose: Convert an instanceMatrix (kxn matrix where n is number of voxels, k is number of instances)
%             to a cell array of k classes of instances
%
function [instances stimulusClasses] = instanceMatrix2instance(instanceMatrix,stimVals)

instances = {};

% check arguments
if ~any(nargin == [2])
  help instanceMatrix2instance
  return
end

% get number of stimvals
k = length(stimVals);

% check against instanceMatrix
kxn = size(instanceMatrix);
% check dimensions
if k ~= kxn(1)
  if k == kxn(2)
    disp(sprintf('(instanceMatrix2instances) instanceMatrix (%ix%i) appears to have instances in cols, taking transpose',kxn(1),kxn(2)));
    instanceMatrix = instanceMatrix';
    kxn = size(instanceMatrix);
  else
    disp(sprintf('(instanceMatrix2instance) instanceMatrix (%ix%i) does not have matching number of instances: %i\n',kxn(1),kxn(2),k));
    return
  end
end

% get n
n = kxn(2);

% unique stimVals
stimulusClasses = unique(stimVals);
nStimulusClasses = length(stimulusClasses);

% make into instances cell array
for iClass = 1:nStimulusClasses
  instances{iClass} = instanceMatrix(find((stimulusClasses(iClass) == stimVals)),:);
end

