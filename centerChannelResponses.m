% centerChannelResponses.m
%
%      usage: channelResponse = centerChannelResponses(channelResponses,stimVals,channel)
%         by: justin gardner
%       date: 09/07/16
%    purpose: centers channelResponse 
%
function [shiftedResponse shiftedStimVals] = centerChannelResponses(channelResponses,stimVals,channel)

% check arguments
if nargin == 0
  help centerChannelResponses
  return
end

% check for stimulus values matching channel preferences
if ~isempty(setdiff(unique(stimVals(:)), channel.channelPref(:)))
  disp(sprintf('(centerChannelOutput) Could not shift channel response for averaging since stimulus values do not match channel preferences'));
  keyboard
  return
end

% which stimulus value to center to
nChannels = channel.info.numFilters;
centerStimVal = floor(nChannels/2);

% get shifted response
shiftedResponse = zeros(size(channelResponses,1),size(channelResponses,2));
for iResponse = 1:size(channelResponses,1)
  % then see how much we have to circularly shift to get the response aligned. 
  shiftVal = centerStimVal-find(stimVals(iResponse) == channel.channelPref) + 1;
  % shift the response
  thisResponse = channelResponses(iResponse,:);
  shiftedResponse(iResponse,:) = circshift(thisResponse(:),shiftVal);
end

if nargout >= 2
  % get the stimulus values (shifted to center)
  shiftedStimVals = circshift(channel.channelPref(:),centerStimVal);
  shiftedStimVals(1:centerStimVal) = shiftedStimVals(1:centerStimVal)-channel.span;
end

