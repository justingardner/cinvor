% channelNoiseModelTest.m
%
%      usage: noiseModelOutput = channelNoiseModelTest(instanceMatrix,stimVals,channel)
%         by: justin gardner
%       date: 09/08/16
%    purpose: makes likelihood functions using noise model created by channelNoiseModelFit
%
function noiseModelOutput = channelNoiseModelTest(instanceMatrix,stimVals,channel)

noisModelOutput = [];
% check arguments
if nargin < 3
  help channelNoiseModelTest
  return
end

nInstances = size(instanceMatrix,1);
% over instances
for iInstance = 1:nInstances
  % now for all parameter values (like orientation)
  for iStim = 1:length(channel.spanValues)
    % compute the expected voxel pattern
    expectedResponse = channel.spanResponse(iStim,:)*channel.channelWeights;
    noiseModelOutput.likelihood(iInstance,iStim) = mvnpdf(instanceMatrix(iInstance,:),expectedResponse,channel.noiseModel.covar);
  end
  % normalize to probability
  noiseModelOutput.likelihood(iInstance,:) = noiseModelOutput.likelihood(iInstance,:)/sum(noiseModelOutput.likelihood(iInstance,:));
end

% now recenter and average
centerIndex = round(length(channel.spanValues)/2);
for iInstance = 1:nInstances
  % find closest value in span for this stimval
  [~,spanIndex] = min(abs(channel.spanValues-stimVals(iInstance)));
  shiftedLikelihood(iInstance,:) = circshift(noiseModelOutput.likelihood(iInstance,:)',centerIndex-spanIndex);
end

% calculate mean and standard error of test set (throwing out any nans that we get)
noiseModelOutput.meanLikelihood = nanmean(shiftedLikelihood);
noiseModelOutput.steLikelihood = nanstd(shiftedLikelihood)./sqrt(sum(~isnan(shiftedLikelihood)));

