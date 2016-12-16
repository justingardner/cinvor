% dispChannelOutput.m
%
%      usage: dispChannelOutput(channelOutput,channel)
%         by: justin gardner
%       date: 09/07/16
%    purpose: Displays channel output returned by testChannels
%
function statsStr = dispChannelOutput(channelOutput,channel,varargin)

% default stats string
statsStr = '';

% check arguments
if nargin < 2
  help dispChannelOutput
  return
end

% get argments
getArgs(varargin,{'MarkerFaceColor=k','MarkerEdgeColor=w','MarkerSize=8','likelihood=0'});

if likelihood
  if ~isfield(channelOutput,'noiseModel')
    disp(sprintf('(dispChannelOutput) No noise model was calculated'));
    return
  end
  % get x-axis centered
  halfSpanLength = round(length(channel.spanValues)/2);
  spanValues = circshift(channel.spanValues(:),halfSpanLength);
  spanValues(1:halfSpanLength) = spanValues(1:halfSpanLength)-channel.span;

  % plot the values
  myerrorbar(spanValues,channelOutput.noiseModel.meanLikelihood,'yError',channelOutput.noiseModel.steLikelihood,'MarkerFaceColor',MarkerFaceColor,'Symbol=-','LineWidth=2','yErrorBarType=fill');hold on
  % stats string to return
  if isfield(channel,'noiseModel')
    statsStr = sprintf('rho=%0.3f sigma=%0.3f tau=%0.3f',channel.noiseModel.rho,channel.noiseModel.sigma,mean(channel.noiseModel.tau));
  end    
  return
end

% shift the response to recenter them
[shiftedResponse shiftedStimVal] = centerChannelResponses(channelOutput.channelResponse,channelOutput.stimVal,channel);

% and take average
averageResponse = mean(shiftedResponse);
steResponse = std(shiftedResponse)/sqrt(channelOutput.n);

% plot the values
myerrorbar(shiftedStimVal,averageResponse,'yError',steResponse,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerSize',MarkerSize);hold on





