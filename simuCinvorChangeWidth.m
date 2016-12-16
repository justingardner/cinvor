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
%%TO SAVE: save('widthChangeData5','allChannelTuning','nll_vals','T','lowKappaVals','noiseStretch')
clear m;
dispFig = 1;
scenario = 1;
% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=22');
N_trials = 22;
% do each scenario
switch scenario
 case 1
  scenarioName = 'Amplitude modulation only';
  testParameter = 'amplitude';
  testParameterValues = [1 0.8 0.2 0.1];
  testParameterValues = [1/1.73, 1];
  
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
T = 1000;
lowKappaVals = [1,1.5,2,2.5,3];% [1,1.5,2,2.5,3]
allChannelTuning = zeros(length(lowKappaVals),length(testParameterValues),M+1);
nll_vals = zeros(length(lowKappaVals),1);
indStart = zeros(1,length(testParameterValues));
indEnd = zeros(1,length(testParameterValues));
for i = 1:length(testParameterValues)
  indStart(i) = 1 + (i-1)*M;
  indEnd(i) = i*M;
end
fitNoise = false;
noiseStretch = 1.75;
  % get instances from model to train and test
for jj = 1:length(lowKappaVals)
  lowKappa = lowKappaVals(jj);
  for trial = 1:T
    permTrainInstances = [];
    permTestInstances = [];
    % run simulation
    for iValue = 1:length(testParameterValues)
      % set the model parameters
      m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue),'kappa',lowKappa,'noise',(6 + noiseStretch*(lowKappaVals(jj) - 2))/5.78);
      if(iValue > 1)
        m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue),'weighting','fixed','myWeights',m(1).ws,'kappa',lowKappa,'noise',(6 + noiseStretch*(lowKappaVals(jj) - 2))/5.78);
      end
      permTrainInstances = [permTrainInstances getCinvorInstances(m(iValue),e)];
      permTestInstances = [permTestInstances getCinvorInstances(m(iValue),e)];
    end
    for iValue = 1:length(testParameterValues)
      trainInstances = permTrainInstances(indStart(iValue):indEnd(iValue));
      %testInstances = getCinvorInstances(m(3-iValue),e);
      testInstances = permTestInstances(indStart(iValue):indEnd(iValue));
      % compute forward model

      channel = buildChannels(trainInstances,e.stimVals,'fitNoise',fitNoise);
      [avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel,'fitNoise',fitNoise);
      %[channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,m(iValue),e);
      % gather outputs
      thisChannel = [avgTestResponse avgTestResponse(1)];
      allChannelTuning(jj,iValue,:)=allChannelTuning(jj,iValue,:) + 1/(T)*reshape(thisChannel,[1,1,M+1]);
    end
  end
  error_mat = ((reshape(allChannelTuning(jj,:,:),size(realChannelTuning)) - realChannelTuning)./errorChannelTuning).^2/2;
  nll_vals(jj) = sum(sum(error_mat(1,1:M)));
end



% display results of simulation
clf;
h=mlrSmartfig(sprintf('simuCinvor%i',scenario));
set(h,'name',sprintf('Scenario: %s',scenarioName));
channelPrefs = [e.stimVals 180];
%allChannelFit(:,1)=allChannelFit(:,1)*m(1).rangeScaleFac;
for i=1:length(testParameterValues)
  subplot(2,1,i);
  xp= channelPrefs(1):1:channelPrefs(end)*m(1).rangeScaleFac;
  %yp=vonMises(allChannelFit(i,1:4),d2r(xp));
  xpOri=xp/m(1).rangeScaleFac;
  legendVect = [];
  legendLabel = {};
  for jj = 1:length(lowKappaVals)
    l(jj) = plot(channelPrefs, reshape((allChannelTuning(jj,i,:)),[1,9]),'Color',[0.2,jj/length(lowKappaVals),1 - jj/length(lowKappaVals)]); hold on;
    legendVect = [legendVect, l(jj)];
    legendLabel{jj} = ['neural tuning fwhm ', num2str(vonMisesFWHM(lowKappaVals(jj)))];
  end
  %plot(xpOri,yp,[getcolor(i),'--']);
  axis([-10,190,0,0.40]);
  if(i == 1)
    hTitle = title('Channel response functions low contrast');
  else
    hTitle = title('Channel response functions high contrast');
  end
  hXLabel = xlabel('Orientation (deg)');
  hYLabel = ylabel('Channel response function');
  hLegend = legend(legendVect,legendLabel);
  set( gca                       , ...
    'FontName'   , 'Helvetica' );
  set([hLegend,hTitle, hXLabel, hYLabel], ...
      'FontName'   , 'AvantGarde');
  set([hLegend,gca]             , ...
      'FontSize'   , 11           );
  set([hXLabel, hYLabel]  , ...
      'FontSize'   , 12          );
  set( hTitle                    , ...
      'FontSize'   , 14          , ...
      'FontWeight' , 'bold'      );

  set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:0.1:1, ...
    'LineWidth'   , 1         );
end



return;







