% buildChannels.m
%      usage: channel = buildChannels(instances,stimValues,varargin)
%         by: justin gardner and taosheng liu
%       date: 07/25/14
%    purpose: Creates channels using the forward model proposed by Brouwer & Heeger (2009).
%             Input1: instances can be those returned from getInstances (see getInstances). It is 
%             a cell array for each category where each cell is a nxk matrix in which n 
%             is the number of repeats and k is the number of dimensions (e.g. voxels).
%             Input2: stimValues is a vector of stimulus value for each instance class
%             If a list of ROIs is passed in as first argument, will build channel
%             for each ROI (note the recursive call at the beginning).
%             Output:This function returns a structure that contains the
%             channel info., the meaning of some important fields are
%             explained below:
%             spanValues: all possible values in the span of
%             the features space (either 180 deg or 360 deg).
%             spanResponse: the channel responses evoked by each stimlus value in the span.
%             channelPref: the preferred stimulus value by each channel.
%             idealStimVals: specify a list of stimulus values to be used
%             for classification, by default set to the stimulus values in the training data.
%             idealStimResponse: the channel response evoked by the idealStimVals.
%             channelResponse: the theoretical channel response evoked by all the training stimuli.
%             channelWeights: the weight of the channel derived from the channelResponse and the actual instances.


function channel = buildChannels(instances,stimVals,varargin)
channel = [];
% check arguments
if any(nargin == [0])
  help buildChannels
  return
end

instanceFieldName=[]; channelFieldName=[]; model=[]; numFilters=[]; exponent=[]; algorithm=[]; dispChannels=[]; verbose=[];
% parse input arguments
[~,~,preprocessArgs] = getArgs(varargin,{'instanceFieldName=instance','channelFieldName=channel','model=sinFilter','numFilters=8','exponent=7','algorithm=pinv','dispChannels=0','verbose=0'});

% see if we are passed in a cell array of rois. If so, then call buildClassifier
% sequentially on each roi and put the output into the field specified by classField
if isfield(instances{1},instanceFieldName) && isfield(instances{1},'name')
  for iROI = 1:length(instances)
    if ~isfield(instances{iROI}.(instanceFieldName),'instances')
      disp(sprintf('(buildChannels) No instances found in %s for %s',instanceFieldName,instances{iROI}.name));
    else
      % put the output into the roi with the field specified by classField
      disppercent(-inf,sprintf('(buildChannels) Building %s channels for ROI %s',algorithm,instances{iROI}.name));
      instances{iROI}.(channelFieldName) = buildChannels(instances{iROI}.(instanceFieldName).instances,stimVals,varargin{:});
      disppercent(inf);
    end
  end
  channel = instances;
  return
end


% preprocess instances
[instances channel] = preprocessInstances(instances,'args',preprocessArgs);
instanceMatrix=[];
stimValVector=[];

if size(instances, 2)~=length(stimVals)
  error('Number of stimulus values much match the number of classes in instances');
end
for istim=1:length(instances)
  stimValVector=[stimValVector, repmat(stimVals(istim),1,size(instances{istim},1))];
  instanceMatrix=[instanceMatrix; instances{istim}];
end

if max(stimVals)-min(stimVals) <=180
  channel.span=180; multiplier=2; %TSL:note the mulitpler trick (original *2 is actually correct)
else
  channel.span=360; multiplier=1;
end
disp(['(buildChannels) Assume feature space spanned by stimuli/channel is ',num2str(channel.span)]);
channel.spanValues=0:1:channel.span-1;
[channel.spanResponse channel.channelPref]=getChannelResponse(channel.spanValues,multiplier,'model',model,'numFilters',numFilters,'exponent',exponent);
if ~isequal(channel.channelPref, stimVals)
  warning('Channels being built have different preferences than the stimulus. The current implementation is likely incorrect under such a setting');
end
channel.idealStimVals=stimVals;
[channel.idealStimResponse temp]=getChannelResponse(stimVals,multiplier,'model',model,'numFilters',numFilters,'exponent',exponent);

[channel.channelResponse temp]=getChannelResponse(stimValVector,multiplier,'model',model,'numFilters',numFilters,'exponent',exponent);
channel.channelWeights=getChannelWeights(channel.channelResponse, instanceMatrix,'algorithm',algorithm);
% channel.channelWeights=channel.channelWeights./repmat(sum(channel.channelWeights,1),size(channel.channelWeights,1),1); % this will normalize the weights, not sure if it's correct 
channel.info.model=model;
channel.info.numFilters=numFilters;
channel.info.exponent=exponent;
channel.info.algorithm=algorithm;


if dispChannels
  smartfig('Channels','reuse'); clf;
  plot(channel.spanValues, channel.spanResponse,'linewidth',2);
  for i=1:length(channel.channelPref)
    thisLegend{i}=strcat('c',num2str(i));
  end
  legend(thisLegend);
  xlabel('Span in degree');
  ylabel('Response (arb unit)');
  title('Response of each channel to all stim values');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getChannelWeights   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channelWeights = getChannelWeights(channelResponse,instanceMatrix,varargin)

getArgs(varargin,{'algorithm=pinv'});

if strcmp(algorithm,'pinv')
  channelWeights = pinv(channelResponse)*instanceMatrix;
else
  disp(sprintf('(buildChannels:getChannelWeights) Unknown algorithm %s.',algorithm));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getChannelResponse   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [channelResponse channelOrientPref] = getChannelResponse(orientationVec,multiplier,varargin)

getArgs(varargin,{'model=sinFilter','numFilters=8','exponent=7'});

if strcmp(model,'sinFilter')
  [channelResponse channelOrientPref] = sinFilter(orientationVec,multiplier,numFilters,exponent);
else
  disp(sprintf('(docinvor) Unknown filter type: %s',model));
end

%%%%%%%%%%%%%%%%%%%
%%   sinFilter   %%
%%%%%%%%%%%%%%%%%%%
function [filterOut filterOrientPref] = sinFilter(orientation,multiplier,numFilters,filterExponent)

numOrientations=length(orientation);
% get filter phases (evenly spaced)
filterPhase = 0:360/numFilters:359;

% get orientation and filter phase in radians
orientation = d2r(orientation(:)*multiplier);
filterPhase = d2r(filterPhase(:));

% handle multiple phases (i.e. multiple filters)
orientation = repmat(orientation,1,numFilters);
filterPhase = repmat(filterPhase',numOrientations,1);

% sinusoid
filterOut = cos(orientation-filterPhase);

% rectify
filterOut = filterOut.*(filterOut>0);

% apply exponent
filterOut = filterOut.^filterExponent;

% return filterOrientPref (which is just the filterPhase in deg divided by 2)
filterOrientPref = r2d(filterPhase(1,:)/multiplier);
