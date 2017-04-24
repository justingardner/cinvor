% simuCinvor.m
%
%      usage: simuTuningChange(scenario)
%         by: dylan cable
%       date: 7/2016
%    purpose: Simulation of the forward encoding model
%             in which model is trained on one tuning 
%             width and tested on another:
%
%             1 = amplitude simulation
%             2 = noise simulation
%             3 = tuning width simulation
%
function FWHM = simuTuningChange(scenario,firstKappa,secondKappa,noise,varargin)

% check arguments
if ~any(nargin == [1 2 3 4])
  help simuTuningChange
  return
end

% get arguments
getArgs(varargin,'dispFig=0');

% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=24');

% do each scenario
switch scenario
 case 1
  scenarioName = 'Tuning width modulation only';
  testParameter = 'kappa';
  testParameterValues = [firstKappa secondKappa];
end

% display which scenario we are running
disp(scenarioName);
m(1) = setCinvorModel(e,testParameter,testParameterValues(1),'noise',noise);
m(2) = setCinvorModel(e,testParameter,testParameterValues(2),'weighting','fixed','myWeights',m(1).ws,'noise',noise);
% run simulation
for iValue = 1:length(testParameterValues)
  % set the model parameters
  % get instances from model to train and test
  trainInstances = getCinvorInstances(m(1),e);
  testInstances = getCinvorInstances(m(iValue),e);
  % compute forward model
  [channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,m(iValue),e);
  % gather outputs
  allChannelTuning(iValue,:)=thisChannel;
  allChannelFit(iValue,:)=thisFit;
  allR2(iValue)=thisR2;
end

% display results of simulation
if(dispFig)
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
end

FWHM = (allChannelFit(:,5));
return;







