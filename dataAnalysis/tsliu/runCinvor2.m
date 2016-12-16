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
nROI = 2;
totalChannelTuning = zeros(nROI,2,M);
posteriors = cell(nROI,nsubj,2);
rho = zeros(nROI,nsubj,2); 
sigma = zeros(nROI,nsubj,2);
tao = zeros(nROI,nsubj,2);
for ROI = 1:nROI
  for subj = 1:nsubj
    load([subjList{subj},'Anal/',analName]);
    s.lvf.roi{ROI}.instance.instances;
    N = s.lvf.roi{ROI}.n;
    e = specifyCinvorExperiment('stimLevel=8','trialPerStim=1');
    trainInstances = s.lvf.roi{ROI}.instance.instances(1:M);
    testInstances = s.lvf.roi{ROI}.instance.instances(M+1:2*M);

    %testInstances = trainInstances;

    % compute forward model
    channel = buildChannels(trainInstances,e.stimVals);
    rho(ROI,subj,1) = channel.rho;
    sigma(ROI,subj,1) = channel.sigma;
    tao(ROI,subj,1) = mean(channel.tao);
    %[channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,e);
    [avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel);
    % gather outputs
    allChannelTuning(1,:)=avgTestResponse;
    posteriors{ROI,subj,1} = posterior;


    trainInstances = s.lvf.roi{ROI}.instance.instances(M+1:2*M);
    testInstances = s.lvf.roi{ROI}.instance.instances(1:M);
    % compute forward model
    channel = buildChannels(trainInstances,e.stimVals);
    rho(ROI,subj,2) = channel.rho;
    sigma(ROI,subj,2) = channel.sigma;
    tao(ROI,subj,2) = mean(channel.tao);
    %[channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,e);
    [avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel);
    % gather outputs
    allChannelTuning(2,:)=avgTestResponse;
    posteriors{ROI,subj,2} = posterior;
    channelPrefs = channel.channelPref;
    totalChannelTuning(ROI,:,:) = totalChannelTuning(ROI,:,:) + reshape(1.0/(nsubj)*allChannelTuning,size(totalChannelTuning(ROI,:,:)));

  end
end
%append first element to end;

channelPrefs = [channelPrefs 180];
totalChannelTuning = cat(3,totalChannelTuning ,totalChannelTuning(:,:,1));
figure;
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
  title('Channel tuning functions train low test high');
  legend([l(1),l(2)],{'Ipsilateral','Contralateral'});

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
  title('Channel tuning functions train high test low');
  legend([l(1),l(2)],{'Ipsilateral','Contralateral'});
end

subplot(5,2,5);
plot(1:2,reshape(rho(1,:,:),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(rho(1,:,:),[nsubj,2]),1),'bo');
title('rho (Ipsilateral)');
axis([0,3,0,0.5]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('train contrast');

subplot(5,2,6);
plot(1:2,reshape(rho(2,:,:),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(rho(2,:,:),[nsubj,2]),1),'bo');
title('rho (Contralateral)');
axis([0,3,0,0.5]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('train contrast');

subplot(5,2,7);
plot(1:2,reshape(sigma(1,:,:),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(sigma(1,:,:),[nsubj,2]),1),'bo');
title('sigma (Ipsilateral)');
axis([0,3,0,1]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('train contrast');

subplot(5,2,8);
plot(1:2,reshape(sigma(2,:,:),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(sigma(2,:,:),[nsubj,2]),1),'bo');
title('sigma (Contralateral)');
axis([0,3,0,1]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('train contrast');

subplot(5,2,9);
plot(1:2,reshape(tao(1,:,:),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(tao(1,:,:),[nsubj,2]),1),'bo');
title('tao (Ipsilateral)');
axis([0,3,0,1.5]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('train contrast');

subplot(5,2,10);
plot(1:2,reshape(tao(2,:,:),[nsubj,2])','ro');
hold on;
plot(1:2,mean(reshape(tao(2,:,:),[nsubj,2]),1),'bo');
title('tao (Contralateral)');
axis([0,3,0,1.5]);
ax = gca;
set(ax,'XTickLabel',{'','low','high',''});
xlabel('train contrast');

roundedStimVals = round(channelPrefs-0.1)+1;
%rotate through to 90;
figure;
for lowHighInd = 1:2
  for thisROI = 1:2
    M = 8;

    allPosterior = zeros(0,8,180);
    for subj = 1:nsubj
      thisPosterior = posteriors(thisROI,subj,lowHighInd);
      trialPerStim = length(thisPosterior{1}.mean)/M;
      allPosterior = cat(1,allPosterior,reshape(thisPosterior{1}.val,[trialPerStim,M,180]));
    end
    for plotIndex = 1:M
      allPosterior(:,plotIndex,:) = circshift(allPosterior(:,plotIndex,:),[0,0,90 - roundedStimVals(plotIndex)]);
    end
    allPosterior(isnan(allPosterior)) = 0;
    subplot(2,2,2*(lowHighInd-1)+thisROI);
    plot(median(reshape(allPosterior,[size(allPosterior,1)*M, 180]),1));
    if(thisROI == 1)
      if(lowHighInd == 1)
        title('Median Posterior Train Low Test High Contrast Ipsilateral');
      else
        title('Median Posterior Train High Test Low Contrast Ipsilateral');
      end
    else
      if(lowHighInd == 1)
        title('Median Posterior Train Low Test High Contrast Contrateral');
      else
        title('Median Posterior Train High Test Low Contrast Contrateral');
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
    thisPosterior = posteriors(thisROI,subj,lowHighInd);
    trialPerStim = length(thisPosterior{1}.mean)/M;
    allPosterior = cat(1,allPosterior,reshape(thisPosterior{1}.val,[trialPerStim,M,180]));
  end
  for plotIndex = 1:M
    allPosterior(:,plotIndex,:) = circshift(allPosterior(:,plotIndex,:),[0,0,90 - roundedStimVals(plotIndex)]);
    allPosterior(isnan(allPosterior)) = 0;
    subplot(M,2,thisROI + 2*(plotIndex-1));
    plot(median(reshape(allPosterior(:,plotIndex,:),[size(allPosterior,1), 180]),1));
  end
end
