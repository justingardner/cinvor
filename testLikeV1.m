% testLikeV1.m
%
%      usage: testLikeV1()
%         by: justin gardner
%       date: 09/12/16
%    purpose: 
%
function retval = testLikeV1()

% check arguments
if ~any(nargin == [0])
  help testLikeV1
  return
end

% nFold cross-val
nFold = 5;

% load data
dataDir = '~/Desktop';
d = load(fullfile(dataDir,'dataV1'));

% convert two instance cell array / stimVals
[instances stimVals] = instanceMatrix2instance(d.instances,d.svec);

% get the crossVal sets
crossVal = getCrossValSets(instances,'nFold',nFold);
mlrSmartfig('testLikeV1fold','reuse');clf;
      
% now train and test over folds
for iFold = 1:nFold
  % comupte train/test of left visual field channels
  [trainInstances testInstances] = getCrossValInstances(instances,crossVal,iFold);
  channelFold(iFold) = buildChannels(trainInstances,stimVals,'fitNoiseModel=1','noiseModelGridSearchOnly=1','noiseModelFitTolerence',0.1,'noiseModelGridSteps=2');
  outputFold(iFold) = testChannels(testInstances,stimVals,channelFold(iFold),'doClassify=0');

  % display fold
  subplot(1,nFold,iFold);
  dispstr = dispChannelOutput(outputFold(iFold),channelFold(iFold),'likelihood=1');
  title(sprintf('Fold %i/%i\n%s',iFold,nFold,dispstr));
  xlabel('Orientation difference from actual (deg)');
  ylabel('Likelihood (p)');
  drawnow;
end

% combine across folds
[output channel] = combineChannelOutputs(outputFold,channelFold,'combineAveraged=0');

% display
mlrSmartfig('testLikeV1','reuse');clf;
dispstr = dispChannelOutput(output,channel,'likelihood=1');
title(sprintf('Average over %i folds\n%s',iFold,nFold,dispstr));
xlabel('Orientation difference from actual (deg)');
ylabel('Likelihood (p)');

keyboard

