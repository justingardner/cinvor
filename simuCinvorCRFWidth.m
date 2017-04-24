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

%TO SAVE: save('widthSearch7','allWidths','kappaVals','noiseVals','T','stickFilter','allR2')  
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
kappaVals = [1,2,4,7,10]; %for left termnial, [1.5,2,2.5].  %for right terminal, [1.5,2,2.5,3,4]
noiseVals = [0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
allWidths = zeros(length(kappaVals),length(noiseVals));
allR2 = zeros(length(kappaVals),length(noiseVals));
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
        [channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,m(iValue),e,'model',stickFilter);
        % gather outputs
        allChannelTuning(iValue,:)=allChannelTuning(iValue,:) + 1/T*thisChannel;
        allChannelFit(iValue,:)=thisFit;
        allR2(kappaInd,noiseInd)=allR2(kappaInd,noiseInd) + 1/T*thisR2;
      end
    end
    myFit = fitTuningWithVM([allChannelTuning; allChannelTuning],channelPrefs,2,0);
    
    allWidths(kappaInd,noiseInd) =allWidths(kappaInd,noiseInd)+ myFit(1,5);
  end
end
color_array = [105,217,255;
140,232,121;
255,226,145;
232,129,121;
172,133,255];
figure;
legend_array = {};
for i = 1:length(kappaVals)
  toPlot = allR2(i,:) > 0;
  plot(allR2(i,toPlot),allWidths(i,toPlot),'Color',color_array(i,:)/256);
  hold on;
  xlabel('r2');
  ylabel('crf fwhm')
  legend_array{i} = ['neural tuning width = ',num2str(vonMisesFWHM(kappaVals(i)))];
end

legend(legend_array);
axis([0,1,20,60]);
drawPublishAxis;








