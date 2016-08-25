%toSave: save('channelTuning7','totalChannelTuning','subjChannelTuning','whiten','currVF','whitenMeanOnly','whitenStdOnly','trainTestSame')

subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
nsubj=length(subjList);
analName='decon1gIns';
clear allChannelTuning;
ROI = 1;
nsubj=6;
M = 8;
plotMeans = false;
if(plotMeans)
  thisMeans = zeros(2,M);
  for i = 1:M
    thisMeans(1,i) = mean(mean(s.lvf.roi{ROI}.instance.instances{i},1));
    thisMeans(2,i) = mean(mean(s.lvf.roi{ROI}.instance.instances{i+M},1));
  end
  figure
  plot(thisMeans(1,:));
  hold
  plot(thisMeans(2,:),'r');
end
N_trials = 27;
N_folds = 5;
fold_endpoints = round(linspace(1,N_trials+1,N_folds+1));
nROI = 2;
totalChannelTuning = zeros(nROI,2,M);
subjChannelTuning = zeros(nROI,nsubj,2,M);
posteriors = cell(nROI,nsubj,N_folds,2);
rho = zeros(nROI,nsubj,N_folds,2); 
sigma = zeros(nROI,nsubj,N_folds,2);
tao = zeros(nROI,nsubj,N_folds,2);
trainTestCorr = zeros(nROI,M,2);
trainSTD = zeros(nROI,M,2);
channelWeightMean = zeros(nROI,2);
channelWeightStd = zeros(nROI,2);
whiten = false;
whitenMeanOnly = false;
whitenStdOnly = false;
fitNoise = false;
shuffleInstances = false;
trainTestSame = false; %if false, (trail low, test high).  If true, (train low, test low)
currVF = 'rvf';
for subj = 1:nsubj
  for ROI = 1:nROI
    load([subjList{subj},'Anal/',analName]);
    thisStim = getfield(s,currVF);
    instanceMatrix=[];
    for istim=1:length(thisStim.roi{ROI}.instance.instances)
      instanceMatrix=[instanceMatrix; thisStim.roi{ROI}.instance.instances{istim}];
    end
    voxelMeans = mean(instanceMatrix);
    voxelStd = std(instanceMatrix);
    N = thisStim.roi{ROI}.n;
    e = specifyCinvorExperiment('stimLevel=8','trialPerStim=1');
    for fold_num = 1:N_folds
      trainInstances = thisStim.roi{ROI}.instance.instances(1:M);
      if(trainTestSame)
        testInstances = thisStim.roi{ROI}.instance.instances(1:M);
      else
        testInstances = thisStim.roi{ROI}.instance.instances(M+1:2*M);
      end
      for i = 1:M
        trainInstances{i} = trainInstances{i}([1:fold_endpoints(fold_num)-1 fold_endpoints(fold_num+1):26],:);
        testInstances{i} = testInstances{i}(fold_endpoints(fold_num):fold_endpoints(fold_num+1)-1,:);
        if(whiten)
          thisLengthTrain = size(trainInstances{i},1);
          thisLengthTest = size(testInstances{i},1);
          trainInstances{i} = (trainInstances{i} - repmat(voxelMeans,thisLengthTrain,1))./repmat(voxelStd,thisLengthTrain,1);
          testInstances{i} = (testInstances{i} - repmat(voxelMeans,thisLengthTest,1))./repmat(voxelStd,thisLengthTest,1);
        elseif(whitenMeanOnly)
          thisLengthTrain = size(trainInstances{i},1);
          thisLengthTest = size(testInstances{i},1);
          trainInstances{i} = (trainInstances{i} - repmat(voxelMeans,thisLengthTrain,1));
          testInstances{i} = (testInstances{i} - repmat(voxelMeans,thisLengthTest,1));
        elseif(whitenStdOnly)
          thisLengthTrain = size(trainInstances{i},1);
          thisLengthTest = size(testInstances{i},1);
          trainInstances{i} = (trainInstances{i} - repmat(voxelMeans,thisLengthTrain,1))./repmat(voxelStd,thisLengthTrain,1) + repmat(voxelMeans,thisLengthTrain,1);
          testInstances{i} = (testInstances{i} - repmat(voxelMeans,thisLengthTest,1))./repmat(voxelStd,thisLengthTest,1) + repmat(voxelMeans,thisLengthTest,1);
        end
        trainTestCorr(ROI,i,1) = trainTestCorr(ROI,i,1) + 1.0/(nsubj*N_folds)*corr(mean(trainInstances{i})',mean(testInstances{i})');
        trainSTD(ROI,i,1) = trainSTD(ROI,i,1) + 1.0/(nsubj*N_folds)*mean(std(trainInstances{i}));
      end
      %shuffle test instances
      if(shuffleInstances)
        instanceMatrix=[];
        for istim=1:length(trainInstances)
          instanceMatrix=[instanceMatrix; testInstances{istim}];
        end
        trialPerStim = size(instanceMatrix,1)/M;
        ordering = randperm(size(instanceMatrix,1));
        instanceMatrix = instanceMatrix(ordering,:);
        for istim=1:length(testInstances)
          startIndex = 1 + (istim-1)*trialPerStim;
          endIndex = (istim)*trialPerStim;
          testInstances{istim} = instanceMatrix(startIndex:endIndex,:);
        end
      end
      %testInstances = trainInstances;

      % compute forward model
      channel = buildChannels(trainInstances,e.stimVals,'fitNoise',fitNoise);
      channelWeightMean(ROI,1) = channelWeightMean(ROI,1) + 1.0/(nsubj*N_folds)*mean(mean(channel.channelWeights));
      channelWeightStd(ROI,1) = channelWeightStd(ROI,1) + 1.0/(nsubj*N_folds)*mean(std(channel.channelWeights));
      if(fitNoise)
        rho(ROI,subj,fold_num,1) = channel.rho;
        sigma(ROI,subj,fold_num,1) = channel.sigma;
        tao(ROI,subj,fold_num,1) = mean(channel.tao);
      end
      %[channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,e);
      [avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel,'fitNoise',fitNoise);
      % gather outputs
      allChannelTuning(1,:)=avgTestResponse;
      posteriors{ROI,subj,fold_num,1} = posterior;


      trainInstances = thisStim.roi{ROI}.instance.instances(M+1:2*M);
      if(trainTestSame)
        testInstances = thisStim.roi{ROI}.instance.instances(M+1:2*M);
      else
        testInstances = thisStim.roi{ROI}.instance.instances(1:M);
      end
      for i = 1:M
        trainInstances{i} = trainInstances{i}([1:fold_endpoints(fold_num)-1 fold_endpoints(fold_num+1):26],:);
        testInstances{i} = testInstances{i}(fold_endpoints(fold_num):fold_endpoints(fold_num+1)-1,:);
        if(whiten)
          thisLengthTrain = size(trainInstances{i},1);
          thisLengthTest = size(testInstances{i},1);
          trainInstances{i} = (trainInstances{i} - repmat(voxelMeans,thisLengthTrain,1))./repmat(voxelStd,thisLengthTrain,1);
          testInstances{i} = (testInstances{i} - repmat(voxelMeans,thisLengthTest,1))./repmat(voxelStd,thisLengthTest,1);
        elseif(whitenMeanOnly)
          thisLengthTrain = size(trainInstances{i},1);
          thisLengthTest = size(testInstances{i},1);
          trainInstances{i} = (trainInstances{i} - repmat(voxelMeans,thisLengthTrain,1));
          testInstances{i} = (testInstances{i} - repmat(voxelMeans,thisLengthTest,1));
        elseif(whitenStdOnly)
          thisLengthTrain = size(trainInstances{i},1);
          thisLengthTest = size(testInstances{i},1);
          trainInstances{i} = (trainInstances{i} - repmat(voxelMeans,thisLengthTrain,1))./repmat(voxelStd,thisLengthTrain,1) + repmat(voxelMeans,thisLengthTrain,1);
          testInstances{i} = (testInstances{i} - repmat(voxelMeans,thisLengthTest,1))./repmat(voxelStd,thisLengthTest,1) + repmat(voxelMeans,thisLengthTest,1);
        end
        trainTestCorr(ROI,i,2) = trainTestCorr(ROI,i,2) + 1.0/(nsubj*N_folds)*corr(mean(trainInstances{i})',mean(testInstances{i})');
        trainSTD(ROI,i,2) = trainSTD(ROI,i,2) + 1.0/(nsubj*N_folds)*mean(std(trainInstances{i}));
      end
      % compute forward model
      channel = buildChannels(trainInstances,e.stimVals,'fitNoise',fitNoise);
      channelWeightMean(ROI,2) = channelWeightMean(ROI,2) + 1.0/(nsubj*N_folds)*mean(mean(channel.channelWeights));
      channelWeightStd(ROI,2) = channelWeightStd(ROI,2) + 1.0/(nsubj*N_folds)*mean(std(channel.channelWeights));
      if(fitNoise)
        rho(ROI,subj,fold_num,2) = channel.rho;
        sigma(ROI,subj,fold_num,2) = channel.sigma;
        tao(ROI,subj,fold_num,2) = mean(channel.tao);
      end
      %[channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,e);
      [avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel,'fitNoise',fitNoise);
      % gather outputs
      allChannelTuning(2,:)=avgTestResponse;
      posteriors{ROI,subj,fold_num,2} = posterior;
      channelPrefs = channel.channelPref;
      totalChannelTuning(ROI,:,:) = totalChannelTuning(ROI,:,:) + reshape(1.0/(nsubj*N_folds)*allChannelTuning,size(totalChannelTuning(ROI,:,:)));
      subjChannelTuning(ROI,subj,:,:) = subjChannelTuning(ROI,subj,:,:) + reshape(1.0/N_folds*allChannelTuning,size(subjChannelTuning(ROI,subj,:,:)));
    end
  end
end
errorChannelTuning = reshape(std(subjChannelTuning,[],2),[nROI,2,M])./sqrt(nsubj);
%append first element to end;
channelPrefs = [channelPrefs 180];
totalChannelTuning = cat(3,totalChannelTuning ,totalChannelTuning(:,:,1));
errorChannelTuning = cat(3,errorChannelTuning ,errorChannelTuning(:,:,1));
channelFunction = [channel.idealStimResponse(5,:) 0];
if(strcmp(currVF,'rvf'))
  save('realDataRight','totalChannelTuning','errorChannelTuning'); 
else
  save('realDataLeft','totalChannelTuning','errorChannelTuning');
end
if(true)
  lowHighInd = 1;
  h=mlrSmartfig(sprintf('simuCinvor'));
  set(h,'name',sprintf('Scenario:'));
  subplot(5,2,1:2);
  rangeScaleFac = 2;
  allChannelFit = fitTuningWithVM(reshape(totalChannelTuning(:,lowHighInd,:),[nROI,size(totalChannelTuning,3)]),channelPrefs,rangeScaleFac,0);
  allChannelFit(:,1)=allChannelFit(:,1)*rangeScaleFac;
  for i=1:2
    %channelPrefs = e.stimVals;
    xp= channelPrefs(1):1:channelPrefs(end)*rangeScaleFac;
    yp=vonMises(allChannelFit(i,1:4),d2r(xp));
    xpOri=xp/rangeScaleFac;
    plot(channelPrefs, reshape(totalChannelTuning(i,lowHighInd,:),[1,size(totalChannelTuning,3)]),[getcolor(i),'o']); hold on
    l(i) = plot(xpOri,yp,[getcolor(i),'--']);
  end
  %l(3) = plot(channelPrefs,channelFunction,'b');
  title('Channel tuning functions low contrast');
  xlabel('angle');
  ylabel('percent of maximum channel activation');
  %legend([l(1),l(2),l(3)],{'Ipsilateral','Contralateral','Channel Function'});

  lowHighInd = 2;
  subplot(5,2,3:4);
  rangeScaleFac = 2;
  allChannelFit = fitTuningWithVM(reshape(totalChannelTuning(:,lowHighInd,:),[nROI,size(totalChannelTuning,3)]),channelPrefs,rangeScaleFac,0);
  allChannelFit(:,1)=allChannelFit(:,1)*rangeScaleFac;
  for i=1:2
    %channelPrefs = e.stimVals;
    xp= channelPrefs(1):1:channelPrefs(end)*rangeScaleFac;
    yp=vonMises(allChannelFit(i,1:4),d2r(xp));
    xpOri=xp/rangeScaleFac;
    plot(channelPrefs, reshape(totalChannelTuning(i,lowHighInd,:),[1,size(totalChannelTuning,3)]),[getcolor(i),'o']); hold on
    l(i) = plot(xpOri,yp,[getcolor(i),'--']);
  end
  %l(3) = plot(channelPrefs,channelFunction,'b');
  title('Channel tuning functions high contrast');
  xlabel('angle');
  ylabel('percent of maximum channel activation');
  %legend([l(1),l(2),l(3)],{'Ipsilateral','Contralateral','Channel Function'});
end

subplot(5,2,5);
plot(1:2,reshape(mean(rho(1,:,:,:),3),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(mean(rho(1,:,:,:),3),[nsubj,2]),1),'bo');
title('rho (Ipsilateral)');
axis([0,3,0,0.5]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('contrast');

subplot(5,2,6);
plot(1:2,reshape(mean(rho(2,:,:,:),3),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(mean(rho(2,:,:,:),3),[nsubj,2]),1),'bo');
title('rho (Contralateral)');
axis([0,3,0,0.5]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('contrast');

subplot(5,2,7);
plot(1:2,reshape(mean(sigma(1,:,:,:),3),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(mean(sigma(1,:,:,:),3),[nsubj,2]),1),'bo');
title('sigma (Ipsilateral)');
axis([0,3,0,1]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('contrast');

subplot(5,2,8);
plot(1:2,reshape(mean(sigma(2,:,:,:),3),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(mean(sigma(2,:,:,:),3),[nsubj,2]),1),'bo');
title('sigma (Contralateral)');
axis([0,3,0,1]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('contrast');

subplot(5,2,9);
plot(1:2,reshape(mean(tao(1,:,:,:),3),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(mean(tao(1,:,:,:),3),[nsubj,2]),1),'bo');
title('tao (Ipsilateral)');
axis([0,3,0,1.5]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('contrast');

subplot(5,2,10);
plot(1:2,reshape(mean(tao(2,:,:,:),3),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(mean(tao(2,:,:,:),3),[nsubj,2]),1),'bo');
title('tao (Contralateral)');
axis([0,3,0,1.5]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('contrast');

if(fitNoise)
  roundedStimVals = round(channelPrefs-0.1)+1;
  %rotate through to 90;
  figure;
  for lowHighInd = 1:2
    for thisROI = 1:2
      M = 8;

      allPosterior = zeros(0,8,180);
      for subj = 1:nsubj
        for fold_num = 1:N_folds
          thisPosterior = posteriors(thisROI,subj,fold_num,lowHighInd);
          trialPerStim = length(thisPosterior{1}.mean)/M;
          allPosterior = cat(1,allPosterior,reshape(thisPosterior{1}.val,[trialPerStim,M,180]));
        end
      end
      for plotIndex = 1:M
        allPosterior(:,plotIndex,:) = circshift(allPosterior(:,plotIndex,:),[0,0,90 - roundedStimVals(plotIndex)]);
      end
      allPosterior(isnan(allPosterior)) = 0;
      subplot(2,2,2*(lowHighInd-1)+thisROI);
      plot(median(reshape(allPosterior,[size(allPosterior,1)*M, 180]),1));
      if(thisROI == 1)
        if(lowHighInd == 1)
          title('Median Posterior Low Contrast Ipsilateral');
        else
          title('Median Posterior High Contrast Ipsilateral');
        end
      else
        if(lowHighInd == 1)
          title('Median Posterior Low Contrast Contrateral');
        else
          title('Median Posterior High Contrast Contrateral');
        end
      end
    end
  end

  figure;
  for thisROI = 1:2
    M = 8;
    lowHighInd = 2;

    allPosterior = zeros(0,8,180);
    for subj = 1:nsubj
      for fold_num = 1:N_folds
        thisPosterior = posteriors(thisROI,subj,fold_num,lowHighInd);
        trialPerStim = length(thisPosterior{1}.mean)/M;
        allPosterior = cat(1,allPosterior,reshape(thisPosterior{1}.val,[trialPerStim,M,180]));
      end
    end
    for plotIndex = 1:M
      allPosterior(:,plotIndex,:) = circshift(allPosterior(:,plotIndex,:),[0,0,90 - roundedStimVals(plotIndex)]);
      allPosterior(isnan(allPosterior)) = 0;
      subplot(M,2,thisROI + 2*(plotIndex-1));
      plot(median(reshape(allPosterior(:,plotIndex,:),[size(allPosterior,1), 180]),1));
    end
  end

  mean_posterior = zeros(nROI,nsubj,2,M);
  std_posterior = zeros(nROI,nsubj,2,M);
  for ROI = 1:nROI
    for subj = 1:nsubj
      for condi = 1:2
        all_mean = zeros(0,M);
        all_std = zeros(0,M);
        for fold_num = 1:N_folds
          thisPosterior = posteriors(ROI,subj,fold_num,condi);
          trialPerStim = length(thisPosterior{1}.mean)/M;
          all_mean = [all_mean; reshape((thisPosterior{1}.mean),[trialPerStim,M])];
          all_std = [all_std; reshape((thisPosterior{1}.std),[trialPerStim,M])];
        end
        mean_posterior(ROI,subj,condi,:) = circ_mean(all_mean/180*(2*pi))/(2*pi)*180;
        std_posterior(ROI,subj,condi,:) = mean(all_std);
      end
    end
  end
end