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
function simuCinvor(scenario,varargin)

% check arguments
if ~any(nargin == [1 2])
  help simuCinvor
  return
end

% get arguments
getArgs(varargin,'dispFig=1');

% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=24');

% do each scenario
switch scenario
 case 1
  scenarioName = 'Amplitude modulation only';
  testParameter = 'amplitude';
  testParameterValues = [1 0.8 0.2 0.1];
  
 case 2
  scenarioName = 'Noise modulation only';
  testParameter = 'noise';
  testParameterValues=[0.09 0.12 0.14 0.18];
  testParameterValues=[0.09 0.12 0.14 0.54];
  testParameterValues=[0.4 0.09];
        
 case 3
  scenarioName = 'Tuning width modulation only';
  testParameter = 'kappa';
  testParameterValues = [4 3 2 1];
  testParameterValues = [1 5];
end

% display which scenario we are running
disp(scenarioName);

% run simulation
for iValue = 1:length(testParameterValues)
  % set the model parameters
  m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue));
  % get instances from model to train and test
  trainInstances = getCinvorInstances(m(iValue),e);
  testInstances = getCinvorInstances(m(iValue),e);
  % compute forward model
  [channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,m(iValue),e);
  % gather outputs
  allChannelTuning(iValue,:)=thisChannel;
  allChannelFit(iValue,:)=thisFit;
  allR2(iValue)=thisR2;
end

% display results of simulation
h=mlrSmartfig(sprintf('simuCinvor%i',scenario));
set(h,'name',sprintf('Scenario: %s',scenarioName));
subplot(5,2,1:2);
allChannelFit(:,1)=allChannelFit(:,1)*m(1).rangeScaleFac;
for i=1:length(testParameterValues)
  xp= channelPrefs(1):1:channelPrefs(end)*m(1).rangeScaleFac;
  yp=vonMises(allChannelFit(i,1:4),d2r(xp));
  xpOri=xp/m(1).rangeScaleFac;
  plot(channelPrefs, allChannelTuning(i,:),'o'); hold on
  plot(xpOri,yp,[getcolor(i),'--']);
end
title('Channel tuning functions');

subplot(5,2,3);
bar(r2d(allChannelFit(:,1)));
set(gca,'xticklabel',cellArray(testParameterValues));
title('mu');

subplot(5,2,4);
bar(allChannelFit(:,3));
set(gca,'xticklabel',cellArray(testParameterValues));
title('baseline');

subplot(5,2,5);
bar(allChannelFit(:,4));
set(gca,'xticklabel',cellArray(testParameterValues));
title('Amplitude');

subplot(5,2,6);
bar(allChannelFit(:,5));
set(gca,'xticklabel',cellArray(testParameterValues));
title('FWHM');

subplot(5,2,7);
bar(allR2);
set(gca,'xticklabel',cellArray(testParameterValues));
title('r2');

subplot(5,2,8);
plot(allR2,allChannelFit(:,3),'o-'); hold on
title('baseline vs. r2')

subplot(5,2,9);
plot(allR2,allChannelFit(:,4),'o-'); hold on
title('amplitude vs. r2');

subplot(5,2,10);
plot(allR2,allChannelFit(:,5),'*--'); 
title('FWHM vs. r2');

return;







