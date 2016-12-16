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
%function simuCinvor(scenario,varargin)

% check arguments
%if ~any(nargin == [1 2])
%  help simuCinvor
%  return
%end

% get arguments
%getArgs(varargin,'dispFig=1');

%TO SAVE: save('boxSearch9','allErrors','kappaVals','noiseVals','T','left','both','trainHighOnly')  
clear m;
dispFig = 1;
scenario = 1;
% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=21');

% do each scenario
switch scenario
 case 1
  scenarioName = 'Amplitude modulation only';
  testParameter = 'amplitude';
  testParameterValues = [1 0.8 0.2 0.1];
  testParameterValues = [1/1.73,1];
  
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
trainHighOnly = true;
left = false;
fitNoise = 0;
both = true;
M = 8;
T = 30;
kappaVals = [1.75,2,2.25,2.5,3,4]; %for left termnial, [1.5,2,2.5].  %for right terminal, [1.5,2,2.5,3,4]
noiseVals = [7.25,7.5,7.75,8,8.25,8.5,8.75,9,9.25,9.5,10];
kappaVals = [4];
noiseVals = [8,10];
allErrors = zeros(length(kappaVals),length(noiseVals));
for kappaInd = 1:length(kappaVals)
  for noiseInd = 1:length(noiseVals)
    allChannelTuning = zeros(length(testParameterValues),M+1);
      
    for trial = 1:T
      % run simulation
      for iValue = 1:length(testParameterValues)
        % set the model parameters
        m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue),'kappa',kappaVals(kappaInd),'noise',noiseVals(noiseInd));
        if(iValue > 1)
          m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue),'weighting','fixed','myWeights',m(1).ws,'kappa',kappaVals(kappaInd),'noise',noiseVals(noiseInd));
        end
      end
      for iValue = 1:length(testParameterValues)
        trainInstances = getCinvorInstances(m(iValue),e);
        %testInstances = getCinvorInstances(m(3-iValue),e);
        testInstances = getCinvorInstances(m(iValue),e);
        % compute forward model
        channel = buildChannels(trainInstances,e.stimVals,'fitNoise',fitNoise);
        [avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel,'fitNoise',fitNoise);
        thisChannel = [avgTestResponse avgTestResponse(1)];
        % gather outputs
        allChannelTuning(iValue,:)=allChannelTuning(iValue,:) + 1/T*thisChannel;
      end
    end
    difference = allChannelTuning-realChannelTuning;
    if(trainHighOnly)
      thisError = (sum(sum(abs(difference(2,1:M)).^2))/(2*M)).^(1/2);
      keyboard
    else
      thisError = (sum(sum(abs(difference(:,1:M)).^2))/(2*M)).^(1/2);
    end
    allErrors(kappaInd,noiseInd) = thisError;
  end
end








