% testChannels.m
%      usage: channelOutput = testChannels(instances,stimValues,channel,varargin)
%         by: justin gardner and taosheng liu
%       date: 07/25/14
%    purpose: Test channel response using the forward model proposed by Brouwer & Heeger (2009).
%
%             instances: instances can be those returned from getInstances (see getInstances), should be
%             the test instance from an independent section of the data (i.e., not the ones
%             used to build the channel).
%             stimValues: stimValues is a vector of stimulus value for each instance class
%             channel: channel is a struct returned by buildChannels
%
%             If a list of ROIs is passed in as first argument, will test channels
%             in each ROI (note the recursive call at the beginning). 
%
%             Outputs are:
%
%             With single return argumnet, a struct channelOutput:
%             
%             channelOutput = testChannels(instances,stimVals,channel);
% 
%             channelOutput.channelResponse has the channel response for each test instance
%             channelOutput.n = number of test responses
%             channelOuptut.r2 = the r2 of the channel model fit
%
%             You can display the channelOutput using:
%
%             dispChannelOutput(channelOutput,channel);
%
%             and comine output for different runs of testChannels using combineChannelOutput.
%
%             With multiple return arguments (old style, returns the following fields):
%               avgTestResponse: the average channel response to the test instances
%               r2: the r2 value for predicting voxel responses, different methods are used
%               classifyCorrTotal: the number of correct and total classifications using the channel response.
%               stimValVector: the actual stimulus value for each test instance
%               predStimVal: the predicted stimulus value for each test instance using the channel response

function [avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(instances,stimVals,channel,varargin)

% classifier = [];
% check arguments
if any(nargin == [0])
  help testChannels
  return
end

% parse input arguments
getArgs(varargin,{'instanceFieldName=instance','channelFieldName=channel','verbose=0','fitNoise=1','doClassify=1'});

if isfield(instances{1},instanceFieldName) && isfield(instances{1},'name')
  for iROI = 1:length(instances)
    if ~isfield(instances{iROI}.(instanceFieldName),'instances')
      disp(sprintf('(testChannels) No instances found in %s for %s',instanceFieldName,instances{iROI}.name));
    elseif ~isfield(channel{iROI},channelFieldName) 
      disp(sprintf('(testChannels) No %s field found in %s. In this case the second argument must be a list of ROIs with precomputed channels.',channelFieldName,instances{iROI}.name));
    else
      [avgTestResponse(iROI,:) r2(iROI) classifyCorrTotal(iROI,:) stimValVector(iROI,:) predStimVal(iROI,:)]= testChannels(instances{iROI}.(instanceFieldName).instances,stimVals,channel{iROI}.(channelFieldName));
    end
  end
  return
end

% preprocess instances, using the setting from the trained channel
instances = preprocessInstances(instances, channel);

% create instance matrix
instanceMatrix=[]; stimValVector=[]; %stimClassVec=[];
for istim=1:length(instances)
  stimValVector=[stimValVector, repmat(stimVals(istim),1,size(instances{istim},1))];
  instanceMatrix=[instanceMatrix; instances{istim}];
end

if(fitNoise)
  [posterior.val posterior.mean posterior.std] = getPosterior(channel,instanceMatrix)
else
  posterior.val = 0;
  posterior.mean = 0;
  posterior.std = 0;
end
% get channel responses
testChannelResponse=instanceMatrix*pinv(channel.channelWeights); 
% and average
avgTestResponse=getAverageChannelResponse(testChannelResponse, stimValVector, channel.channelPref, channel.span/2);


%Prediction/Identification: correlate channel response with the span response and
%find the max correlation, which is the predicted value for the stimulus
corrWithSpanResp=corr(testChannelResponse',channel.spanResponse');
[maxCorr whichStim]=max(corrWithSpanResp,[],2);
predStimVal=channel.spanValues(whichStim);
%now do the same trick but correlate it with ideal reponse so we can classify the stimulus label
corrWithIdealResp=corr(testChannelResponse',channel.idealStimResponse');
[maxCorr whichStim]=max(corrWithIdealResp,[],2);
classifiedStimval=channel.idealStimVals(whichStim);
classifyCorrTotal=[sum(classifiedStimval==stimValVector) length(stimValVector)];

%this part get predicted instance for each stim class 
for i=1:length(stimVals)
  idx=channel.idealStimVals==stimVals(i);
  thisResp=channel.idealStimResponse(idx,:);
%   predInstances(i,:)=thisResp*channel.channelWeights*max(avgTestResponse);
  predInstances(i,:)=thisResp*channel.channelWeights;
end


% package up
if nargout == 1
  channelOutput.n = length(stimValVector);
  channelOutput.stimVal = stimValVector;
  channelOutput.channelResponse = testChannelResponse;
  channelOutput.averageChannelResponse = avgTestResponse;
end

% get likelihood function
if isfield(channel,'noiseModel')
  channelOutput.noiseModel = channelNoiseModelTest(instanceMatrix,stimValVector,channel);
end


if doClassify 
  %Prediction/Identification: correlate channel response with the span response and
  %find the max correlation, which is the predicted value for the stimulus
  corrWithSpanResp=corr(testChannelResponse',channel.spanResponse');
  [maxCorr whichStim]=max(corrWithSpanResp,[],2);
  predStimVal=channel.spanValues(whichStim);
  %now do the same trick but correlate it with ideal reponse so we can classify the stimulus label
  corrWithIdealResp=corr(testChannelResponse',channel.idealStimResponse');
  [maxCorr whichStim]=max(corrWithIdealResp,[],2);
  classifiedStimval=channel.idealStimVals(whichStim);
  classifyCorrTotal=[sum(classifiedStimval(:)==stimValVector(:)) length(stimValVector)];

  %this part get predicted instance for each stim class 
  for i=1:length(stimVals)
    idx=channel.idealStimVals==stimVals(i);
    thisResp=channel.idealStimResponse(idx,:);
%   predInstances(i,:)=thisResp*channel.channelWeights*max(avgTestResponse);
    predInstances(i,:)=thisResp*channel.channelWeights;
  end
  nclass=size(predInstances,1);
  nvox=size(predInstances,2);

  %Now calculate a r2 value for the actual test data
  allTestInstance=[]; allPredInstance=[];
  for i=1:nclass
    thisTestClass=instances{i}; 
    thisPredClass=repmat(predInstances(i,:),size(thisTestClass,1),1);
    allTestInstance=[allTestInstance; thisTestClass];
    allPredInstance=[allPredInstance; thisPredClass];
    thisClassSSR=0; thisClassSS=0;
    for j=1:size(thisTestClass,1) %loop through trials
      thisTestTrial=thisTestClass(j,:);
      sumOfSquaresResidual = sum((thisTestTrial-thisPredClass(j,:)).^2);
      sumOfSquares = sum((thisTestTrial-mean(thisTestTrial)).^2);
      thisR2(j)=1-sumOfSquaresResidual/sumOfSquares;

      thisClassSSR=thisClassSSR+sumOfSquaresResidual;
      thisClassSS=thisClassSS+sumOfSquares;    
    end
    r2.clAvg(i)=mean(thisR2);
    r2.clAcc(i)= 1-thisClassSSR/thisClassSS;
    %doing everything in a single shot
    thisTestClass=thisTestClass(:);
    thisPredClass=thisPredClass(:);
    ssr=sum((thisTestClass-thisPredClass).^2);
    ss=sum((thisTestClass-mean(thisTestClass)).^2);
    r2.clOneshot(i) =1-ssr/ss;
  end
  allTestInstance=allTestInstance(:);
  allPredInstance=allPredInstance(:);
  ssr=sum((allTestInstance-allPredInstance).^2);
  ss=sum((allTestInstance-mean(allTestInstance)).^2);
  r2.overall=1-ssr/ss;

  %do it for each voxel
  for i=1:nvox
    temp=cellfun(@(x) x(:,i),instances,'UniformOutput',false);
    thisTestVox=cell2mat(temp);
    thisPredVox=repmat(predInstances(:,i)', size(thisTestVox,1),1);
    thisVoxSSR=0; thisVoxSS=0;
    for j=1:size(thisTestVox, 1)
      thisTestTrial=thisTestVox(j,:);
      sumOfSquaresResidual = sum((thisTestTrial-thisPredVox(j,:)).^2);
      sumOfSquares = sum((thisTestTrial-mean(thisTestTrial)).^2);
      thisR2Vox(j)=1-sumOfSquaresResidual/sumOfSquares;
      
      thisVoxSSR=thisVoxSSR+sumOfSquaresResidual;
      thisVoxSS=thisVoxSS+sumOfSquares;   
    end
    r2.voxAvg(i)=mean(thisR2Vox);
    r2.voxAcc(i)= 1-thisVoxSSR/thisVoxSS;
    %doing everything in a single shot
    thisTestVox=thisTestVox(:);
    thisPredVox=thisPredVox(:);
    ssr=sum((thisTestVox-thisPredVox).^2);
    ss=sum((thisTestVox-mean(thisTestVox)).^2);
    r2.voxOneshot(i) =1-ssr/ss;   
  end
  % package up into one single return value 
  if nargout == 1
    channelOutput.r2 = r2;
    channelOutput.classifyCorrect = classifyCorrTotal;
    channelOutput.predictedStimVal = predStimVal;
  end
end

% just to return the right argument
if nargout == 1
  avgTestResponse = channelOutput;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getAverageChannelResponse    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avgChannelResponse=getAverageChannelResponse(testChannelResponse,stimValVector,channelPref,channelCenter)

% check whether the stim values match the channel preferences
if ~isequal(unique(stimValVector), channelPref)
  disp(sprintf('(testChannels) The current shift and average scheme is probably incorrect when channel preferences are different from the stimulus values - skipping.'));
  avgChannelResponse = [];
  return
end

centerIdx=find(channelPref==channelCenter);
if isempty(centerIdx)
  disp(['should have a channel centered on ', num2str(channelCenter)]);
  keyboard;
end
allStimResp=[];
for i=1:length(channelPref)
  theseRespIdx=stimValVector==channelPref(i);
  theseStimResp=testChannelResponse(theseRespIdx,:);
  allStimResp=[allStimResp; circshift(theseStimResp,[0 centerIdx-i])];
end
avgChannelResponse=mean(allStimResp,1);

function [posterior posterior_mean posterior_std] = getPosterior(channel,instanceMatrix)
omegaInv = inv(channel.omega);
N_trials = size(instanceMatrix,1);
posterior = zeros(N_trials,channel.span);
posterior_mean = zeros(N_trials,1);
posterior_std = zeros(N_trials,1);
if(channel.span == 180)
  multiplier = 2;
else
  multiplier = 1;
end
angles = deg2rad(1:multiplier:360);
for i = 1:N_trials
  for j = 1:channel.span
    thisError = instanceMatrix(i,:)' - channel.channelWeights'*channel.spanResponse(j,:)';
    posterior(i,j) = exp(-0.5*thisError'*omegaInv*thisError);
  end
  posterior(i,:) = posterior(i,:)/sum(posterior(i,:)); %normalization
  posterior_mean(i) = circ_mean(angles',posterior(i,:)')/(2*pi)*(360/multiplier);
  posterior_std(i) = circ_std(angles',posterior(i,:)')/(2*pi)*(360/multiplier);
end



