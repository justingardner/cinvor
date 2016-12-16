% simuCinvor.m
%
%      usage: simuCinvor(scenario)
%         by: taosheng liu
%       date: 11/2015
%    purpose: Simulation of the forward encoding model
%             scenario can be:
%
%             1 = amplitude simulation
%             2 = noise simulation
%             3 = tuning width simulation
%
function squared_error = simuCinvorSearch(noise,kappa,alphaParam,realChannelTuning)

% check arguments
%if ~any(nargin == [1 2])
%  help simuCinvor
%  return
%end

% get arguments
%getArgs(varargin,'dispFig=1');
clear m;
dispFig = 1;
scenario = 1;
% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=42');
N_trials = 27;
N_voxels = 78;
% do each scenario
switch scenario
 case 1
  scenarioName = 'Amplitude modulation only';
  testParameter = 'amplitude';
  testParameterValues = [1 0.8 0.2 0.1];
  testParameterValues = [1/1.74 1];
  
 case 2
  scenarioName = 'Noise modulation only';
  testParameter = 'noise';
  testParameterValues=[0.09 0.12 0.14 0.18];
  testParameterValues=[0.09 0.12 0.14 0.54];
  testParameterValues=[0.15 0.1];
        
 case 3
  scenarioName = 'Tuning width modulation only';
  testParameter = 'kappa';
  testParameterValues = [4 3 2 1];
  testParameterValues = [1 4];
end

% display which scenario we are running
disp(scenarioName);
M = 8;
T = 20;
N_folds = 2;
fold_endpoints = round(linspace(1,N_trials+1,N_folds+1));
allChannelTuning = zeros(length(testParameterValues),M+1);
trainTestCorr = zeros(1,M,length(testParameterValues));
trainSTD = zeros(1,M,length(testParameterValues));
channelWeightMean = zeros(1,length(testParameterValues));
channelWeightStd = zeros(1,length(testParameterValues));
currTrialTuning = zeros(T,length(testParameterValues),M+1);
fitNoise = false;
indStart = zeros(1,length(testParameterValues));
indEnd = zeros(1,length(testParameterValues));
whiten = true;
ipsilateral = false;
for i = 1:length(testParameterValues)
	indStart(i) = 1 + (i-1)*M;
	indEnd(i) = i*M;
end
  % get instances from model to train and test
for trial = 1:T
  % run simulation
  permTrainInstances = [];
  for iValue = 1:length(testParameterValues)
    % set the model parameters
    m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue),'noise',noise,'kappa',kappa,'alphaParam',alphaParam);
    if(iValue > 1)
      m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue),'weighting','fixed','myWeights',m(1).ws,'noise',noise,'kappa',kappa,'alphaParam',alphaParam);
    end
    permTrainInstances = [permTrainInstances getCinvorInstances(m(iValue),e)];
  end
  if(ipsilateral)
	  for i = 1:length(permTrainInstances)
	  	permTrainInstances{i} = randn(N_trials,N_voxels);
	  end
	end
  instanceMatrix=[];
  for istim=1:length(permTrainInstances)
    instanceMatrix=[instanceMatrix; permTrainInstances{istim}];
  end
  voxelMeans = mean(instanceMatrix);
  voxelStd = std(instanceMatrix);
  for iValue = 1:length(testParameterValues)
    for fold_num = 1:N_folds
      trainInstances = permTrainInstances(indStart(iValue):indEnd(iValue));
      %testInstances = getCinvorInstances(m(3-iValue),e);
      testInstances = permTrainInstances(indStart(iValue):indEnd(iValue));
      for i = 1:M
        trainInstances{i} = trainInstances{i}([1:fold_endpoints(fold_num)-1 fold_endpoints(fold_num+1):26],:);
        testInstances{i} = testInstances{i}(fold_endpoints(fold_num):fold_endpoints(fold_num+1)-1,:);
        if(whiten)
          thisLengthTrain = size(trainInstances{i},1);
          thisLengthTest = size(testInstances{i},1);
          trainInstances{i} = (trainInstances{i} - repmat(voxelMeans,thisLengthTrain,1))./repmat(voxelStd,thisLengthTrain,1);
          testInstances{i} = (testInstances{i} - repmat(voxelMeans,thisLengthTest,1))./repmat(voxelStd,thisLengthTest,1);
        end
      end
      % compute forward model
      for i = 1:M
        trainTestCorr(1,i,iValue) = trainTestCorr(1,i,iValue) + 1.0/(T)*corr(mean(trainInstances{i})',mean(testInstances{i})');
        trainSTD(1,i,iValue) = trainSTD(1,i,iValue) + 1.0/(T)*mean(std(trainInstances{i}));
      end
      channel = buildChannels(trainInstances,e.stimVals,'fitNoise',fitNoise);
      channelWeightMean(1,iValue) = channelWeightMean(1,iValue) + 1.0/(T)*mean(mean(channel.channelWeights));
      channelWeightStd(1,iValue) = channelWeightStd(1,iValue) + 1.0/(T)*mean(std(channel.channelWeights));
      [avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel,'fitNoise',fitNoise);
      %[channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,m(iValue),e);
      % gather outputs
      thisChannel = [avgTestResponse avgTestResponse(1)];
      allChannelTuning(iValue,:)=allChannelTuning(iValue,:) + 1/(T*N_folds)*thisChannel;
      currTrialTuning(trial,iValue,:)=currTrialTuning(trial,iValue,:) + reshape(1/(N_folds)*thisChannel,[1,1,M+1]);
      %allChannelFit(iValue,:)=thisFit;
      %allR2(iValue)=thisR2;
    end
  end
end
noise = noise
alphaParam = alphaParam
kappa = kappa
squared_error = sum(sum(abs(realChannelTuning - allChannelTuning).^2))
return;


