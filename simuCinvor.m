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
% save('simuCinvorData3','allChannelTuning','currTrialTuning','std_mat','testParameterValues','T','N_folds'); 
% get arguments
%getArgs(varargin,'dispFig=1');
clear m;
dispFig = 1;
scenario = 1;
% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=42');
N_trials = 42;
N_voxels = 78;
% do each scenario
switch scenario
 case 1
  scenarioName = 'Amplitude modulation only';
  testParameter = 'amplitude';
  testParameterValues = [1 0.8 0.2 0.1];
  testParameterValues = [1/1.73 1];
  
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
whiten = false;
ipsilateral = false;
crossTrain = false;
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
    m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue));
    if(iValue > 1)
      m(iValue) = setCinvorModel(e,testParameter,testParameterValues(iValue),'weighting','fixed','myWeights',m(1).ws);
    end
    permTrainInstances = [permTrainInstances getCinvorInstances(m(iValue),e)];
  end
  if(crossTrain)
  	permTestInstances = [];
  	for iValue = 1:length(testParameterValues)
	    permTestInstances = [permTestInstances getCinvorInstances(m(3-iValue),e)];
	  end
  else
  	permTestInstances = [];
  	for iValue = 1:length(testParameterValues)
	    permTestInstances = [permTestInstances getCinvorInstances(m(iValue),e)];
	  end
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
      testInstances = permTestInstances(indStart(iValue):indEnd(iValue));
      for i = 1:M
        trainInstances{i} = trainInstances{i}([1:fold_endpoints(fold_num)-1 fold_endpoints(fold_num+1):N_trials],:);
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
nsubj = 12;
std_mat = reshape(std(currTrialTuning,[],1),size(allChannelTuning))/sqrt(nsubj);
z_score = (realChannelTuning - allChannelTuning)./(std_mat);
squared_error = sum(sum(abs(realChannelTuning - allChannelTuning).^2));
trial_error = sort(sum(sum((currTrialTuning - repmat(reshape(allChannelTuning,[1,length(testParameterValues),M+1]),[T,1,1])).^2,3),2));
num_larger = sum(trial_error > squared_error);

nsubj = 12;
S = 300;
sample_larger = 0;
sample_error = zeros(S,1);
for sample = 1:S
  sampleChannelTuning = zeros(length(testParameterValues),M+1);
  for subj = 1:nsubj
    sampleChannelTuning = sampleChannelTuning + 1.0/nsubj*reshape(currTrialTuning(randi([1,T]),:,:),size(sampleChannelTuning));
  end
  sample_squared_error = sum(sum(abs(allChannelTuning - sampleChannelTuning).^2));
  if(sample_squared_error > squared_error)
    sample_larger = sample_larger + 1;
  end
  sample_error(sample) = sample_squared_error;
end
p_val = sample_larger./S;
%difference = allChannelTuning-realChannelTuning;
%thisError = (sum(sum(abs(difference(:,1:M)).^1.5))/(2*M)).^(1/1.5);


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
  l(1) = plot(channelPrefs, 1/2*(allChannelTuning(i,:) + fliplr(allChannelTuning(i,:))),['k']); hold on;
  errorbar(channelPrefs, 1/2*(allChannelTuning(i,:) + fliplr(allChannelTuning(i,:))),1.96*std_mat(i,:),['ko']);
  %plot(xpOri,yp,[getcolor(i),'--']);
  axis([-10,190,0,0.35]);
  hold on;
  l(2) = plot(channelPrefs,realChannelTuning2(i,:),'or','MarkerFaceColor','r','MarkerEdgeColor',[1,1,1],'MarkerSize',10);
  errbar(channelPrefs,realChannelTuning2(i,:),errorChannelTuning(i,:),'-r');
  drawPublishAxis;
  if(i == 1)
    hTitle = title('Channel response functions low contrast');
  else
    hTitle = title('Channel response functions high contrast');
  end
  hXLabel = xlabel('Orientation (deg)');
  hYLabel = ylabel('Channel response function');
  hLegend = legend([l(1),l(2)],{'Model','Data'});
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
vm_fit_data = fitTuningWithVM(realChannelTuning2,channelPrefs,2,1);

nsubj = 12;
S = 300;
vm_fit_sample = zeros(S,2,5);
for sample = 1:S
  sampleChannelTuning = zeros(length(testParameterValues),M+1);
  for subj = 1:nsubj
    sampleChannelTuning = sampleChannelTuning + 1.0/nsubj*reshape(currTrialTuning(randi([1,T]),:,:),size(sampleChannelTuning));
  end
  vm_fit_sample(sample,:,:) = reshape(fitTuningWithVM(sampleChannelTuning,channelPrefs,2,0),[2,5]);
  %reshape(fitTuningWithVM(sampleChannelTuning,channelPrefs,2,0),[2,5]);
end
vm_fit_mean = reshape(mean(vm_fit_sample,1),[2,5]);
vm_fit_std = reshape(std(vm_fit_sample,1),[2,5]);

figure;
subplot(3,2,1);
plot(1:2,vm_fit_data(:,1),'ro');
hold on;
plot(1:2,vm_fit_mean(:,1),'ko');
legend('data','model');
set(gca,'xtick',[1 2]);
set(gca,'XTickLabel',str2mat('low','high'));
xlabel('contrast');
errorbar(1:2,vm_fit_mean(:,1),1.96*vm_fit_std(:,1),'k');
title('mu');

subplot(3,2,2);
plot(1:2,vm_fit_data(:,2),'ro');
hold on;
plot(1:2,vm_fit_mean(:,2),'ko');
legend('data','model');
set(gca,'xtick',[1 2]);
set(gca,'XTickLabel',str2mat('low','high'));
xlabel('contrast');
errorbar(1:2,vm_fit_mean(:,2),1.96*vm_fit_std(:,2),'k');
title('kappa');

subplot(3,2,3);
plot(1:2,vm_fit_data(:,3),'ro');
hold on;
plot(1:2,vm_fit_mean(:,3),'ko');
legend('data','model');
set(gca,'xtick',[1 2]);
set(gca,'XTickLabel',str2mat('low','high'));
xlabel('contrast');
errorbar(1:2,vm_fit_mean(:,3),1.96*vm_fit_std(:,3),'k');
title('base');

subplot(3,2,4);
plot(1:2,vm_fit_data(:,4),'ro');
hold on;
plot(1:2,vm_fit_mean(:,4),'ko');
legend('data','model');
set(gca,'xtick',[1 2]);
set(gca,'XTickLabel',str2mat('low','high'));
xlabel('contrast');
errorbar(1:2,vm_fit_mean(:,4),1.96*vm_fit_std(:,4),'k');
title('amp');

subplot(3,2,5);
plot(1:2,vm_fit_data(:,5),'ro');
hold on;
plot(1:2,vm_fit_mean(:,5),'ko');
legend('data','model');
set(gca,'xtick',[1 2]);
set(gca,'XTickLabel',str2mat('low','high'));
xlabel('contrast');
errorbar(1:2,vm_fit_mean(:,5),1.96*vm_fit_std(:,5),'k');
title('fwhm');


%allChannelFit(:,1) = [0; 0];
%subplot(5,2,3);
%bar(r2d(allChannelFit(:,1)));
%set(gca,'xticklabel',cellArray(testParameterValues));
%title('mu');

%subplot(5,2,4);
%bar(allChannelFit(:,3));
%set(gca,'xticklabel',cellArray(testParameterValues));
%title('baseline');

%subplot(5,2,5);
%bar(allChannelFit(:,4));
%set(gca,'xticklabel',cellArray(testParameterValues));
%title('Amplitude');

%subplot(5,2,6);
%bar(allChannelFit(:,5));
%set(gca,'xticklabel',cellArray(testParameterValues));
%title('FWHM');

%subplot(5,2,7);
%bar(allR2);
%set(gca,'xticklabel',cellArray(testParameterValues));
%title('r2');

%subplot(5,2,8);
%plot(allR2,allChannelFit(:,3),'o-'); hold on
%title('baseline vs. r2')

%subplot(5,2,9);
%plot(allR2,allChannelFit(:,4),'o-'); hold on
%title('amplitude vs. r2');

%subplot(5,2,10);
%plot(allR2,allChannelFit(:,5),'*--'); 
%title('FWHM vs. r2');

return;







