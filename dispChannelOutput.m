% dispChannelOutput.m
%
%      usage: dispChannelOutput(channelOutput,channel)
%         by: justin gardner
%       date: 09/07/16
%    purpose: Displays channel output returned by testChannels
%
function retval = dispChannelOutput(channelOutput,channel,varargin)

% check arguments
if nargin == [0]
  help dispChannels
  return
end

% get argments
getArgs(varargin,{'MarkerFaceColor=k','MarkerEdgeColor=w','MarkerSize=8'});

% shift the response to recenter them
[shiftedResponse shiftedStimVal] = centerChannelResponses(channelOutput.channelResponse,channelOutput.stimVal,channel);

% and take average
averageResponse = mean(shiftedResponse);
steResponse = std(shiftedResponse)/sqrt(channelOutput.n);

% plot the values
myerrorbar(shiftedStimVal,averageResponse,'yError',steResponse,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerSize',MarkerSize);hold on





