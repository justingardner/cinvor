% combineChannelOutputs.m
%
%      usage: channelOutput = combineChannelOutputs(channelOutputs,channel)
%         by: justin gardner
%       date: 09/07/16
%    purpose: Combine multiple channelOutputs returned from testChannels (for example if you want to average over folds)
%
function combined = combineChannelOutputs(channelOutputs,channel,varargin)

% check arguments
if nargin == 0
  help combineChannelOutputs
  return
end

% make sure we have an array
channelOutputs = channelOutputs(:);

% get arguments
getArgs(varargin,{'combineAveraged=0'});

% initialized combined structure
combined.n = 0;
combined.stimVal = [];
combined.channelResponse = [];

% number of channel outputs to average over
nOutputs = length(channelOutputs);

% average over all instances if called for 
% this shifts so that the channelOutputs are centered to the stimulus value for that trial
if combineAveraged
  for iOutput = 1:nOutputs
    % we will end up with one averaged trial
    channelOutputs(iOutput).n = 1;
    % shift the responses
    shiftedResponse = centerChannelResponses(channelOutputs(iOutput).channelResponse,channelOutputs(iOutput).stimVal,channel);
    % and keep the mean
    channelOutputs(iOutput).channelResponse = mean(shiftedResponse);
    % reset the stimulus value
    channelOutputs(iOutput).stimVal = 0;
  end
end

% now combine channel Outputs
for iOutput = 1:nOutputs
  % number of trials
  combined.n = combined.n + channelOutputs(iOutput).n;
  % stimulus values
  combined.stimVal = [combined.stimVal(:) ; channelOutputs(iOutput).stimVal(:)];
  % channel responses
  combined.channelResponse = [combined.channelResponse ; channelOutputs(iOutput).channelResponse];
  % get the r2 value
  combined.r2.overall(iOutput) = channelOutputs(iOutput).r2.overall;
end

% average the r2 value
combined.r2.steOverall = std(combined.r2.overall)/sqrt(nOutputs);
combined.r2.overall = mean(combined.r2.overall);
