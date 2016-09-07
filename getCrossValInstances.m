% getCrossValInstances.m
%
%      usage: [trainInstances testInstances] = getCrossValInstances(instances,crossVal,iFold)
%         by: justin gardner
%       date: 09/07/16
%    purpose: Returns train and test instances for iFold using structure returned by getCrossValSet
%
function [trainInstances testInstances] = getCrossValInstances(instances,crossVal,iFold)

% check arguments
if ~any(nargin == [3])
  help getCrossValInstances
  return
end

% init values
trainInstances = [];
testInstances = [];

% check if iFold is in range
if (iFold < 1) || (iFold > crossVal.nFold)
  disp(sprintf('(getCrossValInstnaces) Fold %i is out of range [1:%i]\n',iFold,crossVal.nFold));
  return
end

% get train and test instances based on crossVal structure
for iClass = 1:size(instances,2)
  trainInstances{iClass} = instances{iClass}(crossVal.trainInstances{iClass,iFold},:);
  testInstances{iClass} = instances{iClass}(crossVal.testInstances{iClass,iFold},:);
end





