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

%TO SAVE: save('channelWidthSearch2','allWidths','expVals','noiseVals','T','stickFilter','allR2')  
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
  testParameterValues = [1];
  
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
stickFilter= 'sinFilter';
trainHighOnly = true;
left = false;
both = true;
M = 8;
T = 1000;
expVals = [1,3,5,7,15]; %for left termnial, [1.5,2,2.5].  %for right terminal, [1.5,2,2.5,3,4]
noiseVals = [0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
allWidths = zeros(length(expVals),length(noiseVals));
allR2 = zeros(length(expVals),length(noiseVals));
for expInd = 1:length(expVals)
  for noiseInd = 1:length(noiseVals)
    allChannelTuning = zeros(length(testParameterValues),M+1);
      
    for trial = 1:T
      % run simulation
      for iValue = 1:length(testParameterValues)
        % set the model parameters
        m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue),'kappa',4,'noise',noiseVals(noiseInd));
      end
      for iValue = 1:length(testParameterValues)
        trainInstances = getCinvorInstances(m(iValue),e);
        %testInstances = getCinvorInstances(m(3-iValue),e);
        testInstances = getCinvorInstances(m(iValue),e);
        % compute forward model
        [channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,m(iValue),e,'model',stickFilter,'exponent',expVals(expInd));
        % gather outputs
        allChannelTuning(iValue,:)=allChannelTuning(iValue,:) + 1/T*thisChannel;
        allChannelFit(iValue,:)=thisFit;
        allR2(expInd,noiseInd)=allR2(expInd,noiseInd) + 1/T*thisR2;
      end
    end
    myFit = fitTuningWithVM([allChannelTuning; allChannelTuning],channelPrefs,2,0);
    
    allWidths(expInd,noiseInd) =allWidths(expInd,noiseInd)+ myFit(1,5);
  end
end
figure;
for i = 1:length(expVals)
  plot(allR2(i,:),allWidths(i,:),getcolor(i));
  hold on;
  xlabel('r2');
  ylabel('crf fwhm')
end
legend('k=1','k=2','k=4','k=7','k=10');








