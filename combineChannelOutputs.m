% combineChannelOutputs.m
%
%      usage: channelOutput = combineChannelOutputs(channelOutputs,channel)
%         by: justin gardner
%       date: 09/07/16
%    purpose: Combine multiple channelOutputs returned from testChannels (for example if you want to average over folds)
%             Optional argument combineAveraged defaults to true, which means that the channelOutput will
%             first be averaged across all test instances (by appropriately shifting the response based
%             on the stimulus value to center the response, see centerChannelResponses). That way, the
%             final combined channel outputs will have one channel response for each channelOutput passed
%             in and when displayed with dispChannelOutput will calculate error bars over the number
%             of channelOutputs passed in here (that is, if you combine over n subjects you will get error
%             bars over n subjects)
%
function [combinedOutput combinedChannel] = combineChannelOutputs(channelOutputs,channel,varargin)

% check arguments
if nargin == 0
  help combineChannelOutputs
  return
end

% make sure we have an array
channelOutputs = channelOutputs(:);

% get arguments
getArgs(varargin,{'combineAveraged=1'});

% initialized combined structure
combinedOutput.n = 0;
combinedOutput.stimVal = [];
combinedOutput.channelResponse = [];

% number of channel outputs to average over
nOutputs = length(channelOutputs);

% degenerate case - just return
if nOutputs == 1
  combinedOutput = channelOutputs;
end

% average over all instances if called for 
% this shifts so that the channelOutputs are centered to the stimulus value for that trial
if combineAveraged
  for iOutput = 1:nOutputs
    % we will end up with one averaged trial
    channelOutputs(iOutput).n = 1;
    % get which channel to use
    if length(channel) < iOutput
      thisChannel = channel(end);
    else
      thisChannel = channel(iOutput);
    end
    % shift the responses
    shiftedResponse= centerChannelResponses(channelOutputs(iOutput).channelResponse,channelOutputs(iOutput).stimVal,thisChannel);
    % and keep the mean
    channelOutputs(iOutput).channelResponse = mean(shiftedResponse);
    % reset the stimulus value
    channelOutputs(iOutput).stimVal = thisChannel.channelPref(round(thisChannel.info.numFilters/2)+1);
  end
end

% now combine channel Outputs
for iOutput = 1:nOutputs
  % number of trials
  combinedOutput.n = combinedOutput.n + channelOutputs(iOutput).n;
  % stimulus values
  combinedOutput.stimVal = [combinedOutput.stimVal(:) ; channelOutputs(iOutput).stimVal(:)];
  % channel responses
  combinedOutput.channelResponse = [combinedOutput.channelResponse ; channelOutputs(iOutput).channelResponse];
  % get the r2 value
  combinedOutput.r2.overall(iOutput) = channelOutputs(iOutput).r2.overall;
end

% average the r2 value
combinedOutput.r2.steOverall = std(combinedOutput.r2.overall)/sqrt(nOutputs);
combinedOutput.r2.overall = mean(combinedOutput.r2.overall);

% combine noiseModelOutputs
if isfield(channelOutputs,'noiseModel')
  for iOutput = 1:nOutputs
    combinedOutput.noiseModel.likelihood(iOutput,:) = channelOutputs(iOutput).noiseModel.meanLikelihood;
  end
  combinedOutput.noiseModel.meanLikelihood = nanmean(combinedOutput.noiseModel.likelihood);
  combinedOutput.noiseModel.steLikelihood = nanstd(combinedOutput.noiseModel.likelihood)/sqrt(sum(~isnan(combinedOutput.noiseModel.likelihood)));
end

% combine channels
if nargout == 2
  combinedChannel = channel(1);
  if isfield(combinedChannel,'noiseModel')
    combinedChannel.noiseModel = [];
    for iChannel = 1:length(channel)
      % average parameters 
      combinedChannel.noiseModel.rho(iChannel) = channel(iChannel).noiseModel.rho;
      combinedChannel.noiseModel.sigma(iChannel) = channel(iChannel).noiseModel.sigma;
      combinedChannel.noiseModel.tau(iChannel) = mean(channel(iChannel).noiseModel.tau);
    end
    % get the average parameter estimates over all folds
    combinedChannel.noiseModel.rho = mean(combinedChannel.noiseModel.rho);
    combinedChannel.noiseModel.sigma = mean(combinedChannel.noiseModel.sigma);
    combinedChannel.noiseModel.tau = mean(combinedChannel.noiseModel.tau);
  end
end