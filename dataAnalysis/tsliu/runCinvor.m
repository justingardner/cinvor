subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
nsubj=length(subjList);
analName='decon1gIns';
totalChannelTuning = zeros(2,9);
nROI = 1;
nStimList = 1;
stimList = {};
stimList{1} = 'lvf';
stimList{2} = 'rvf';
nsubj=1;
M = 8;
plotMeans = false;
if(plotMeans)
  thisMeans = zeros(2,M);
  for i = 1:M
    thisMeans(1,i) = mean(mean(s.lvf.roi{ROI}.instance.instances{i},1));
    thisMeans(2,i) = mean(mean(s.lvf.roi{ROI}.instance.instances{i+M},1));
  end
  figure
  plot(thisMeans(1,:)*4);
  hold
  plot(thisMeans(2,:),'r');
end
N_trials = 216;
trialPerStim = 27;
rho = zeros(nStimList,nROI,2,nsubj); %stimList/ROI/low, high/subj
sigma = zeros(nStimList,nROI,2,nsubj);
tao = zeros(nStimList,nROI,2,nsubj);
posterior_mean = zeros(nStimList,nROI,2,nsubj,N_trials);
posterior_std = zeros(nStimList,nROI,2,nsubj,N_trials);
posterior = zeros(nStimList,nROI,2,nsubj,N_trials,180);
totalChannelTuning = zeros(nStimList,nROI,2,M);
for stimNum = 1:nStimList
	for ROI = 1:nROI
		for subj = 1:nsubj
		  load([subjList{subj},'Anal/',analName]);
		  s.lvf.roi{ROI}.instance.instances;
		  N = s.lvf.roi{ROI}.n;
		  e = specifyCinvorExperiment('stimLevel=8','trialPerStim=1');

			trainInstances = s.lvf.roi{ROI}.instance.instances(1:M);
			%shuffle trainInstances
			instanceMatrix=[];
			for istim=1:length(trainInstances)
			  instanceMatrix=[instanceMatrix; trainInstances{istim}];
			end
			ordering = randperm(N_trials);
			instanceMatrix = instanceMatrix(ordering,:);
			for istim=1:length(trainInstances)
				startIndex = 1 + (istim-1)*trialPerStim;
				endIndex = (istim)*trialPerStim;
			  trainInstances{istim} = instanceMatrix(startIndex:endIndex,:);
			end
			testInstances = s.lvf.roi{ROI}.instance.instances(M+1:2*M);
			%testInstances = trainInstances;

			% compute forward model
			channel = buildChannels(trainInstances,e.stimVals);
			rho(stimNum,ROI,1,subj) = channel.rho;
			sigma(stimNum,ROI,1,subj) = channel.sigma;
			tao(stimNum,ROI,1,subj) = mean(channel.tao);
			posterior_mean(stimNum,ROI,1,subj,:) = channel.posterior_mean;
			posterior_std(stimNum,ROI,1,subj,:) = channel.posterior_std;
			posterior(stimNum,ROI,1,subj,:,:) = channel.posterior;
			%[channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,e);
			[avgTestResponse r2 classifyCorrTotal stimValVector predStimVal] = testChannels(testInstances,e.stimVals,channel);
			% gather outputs
			allChannelTuning(1,:)=avgTestResponse;


			trainInstances = s.lvf.roi{ROI}.instance.instances(M+1:2*M);
			testInstances = s.lvf.roi{ROI}.instance.instances(1:M);

			% compute forward model
			channel = buildChannels(trainInstances,e.stimVals);
			rho(stimNum,ROI,2,subj) = channel.rho;
			sigma(stimNum,ROI,2,subj) = channel.sigma;
			tao(stimNum,ROI,2,subj) = mean(channel.tao);
			posterior_mean(stimNum,ROI,2,subj,:) = channel.posterior_mean;
			posterior_std(stimNum,ROI,2,subj,:) = channel.posterior_std;
			posterior(stimNum,ROI,2,subj,:,:) = channel.posterior;
			%[channelPrefs thisChannel thisFit thisR2] = fitCinvorForwardModel(trainInstances,e.stimVals,testInstances,e.stimVals,e);
			[avgTestResponse r2 classifyCorrTotal stimValVector predStimVal] = testChannels(testInstances,e.stimVals,channel);
			% gather outputs
			allChannelTuning(2,:)=avgTestResponse;

			totalChannelTuning(stimNum,ROI,:,:) = totalChannelTuning(stimNum,ROI,:,:) + reshape((1.0/(nsubj)*allChannelTuning),size(totalChannelTuning(stimNum,ROI,:,:)));

		end
	end
end
if(true)
  h=mlrSmartfig(sprintf('simuCinvor'));
  set(h,'name',sprintf('Scenario:'));
  subplot(5,2,1:2);
  allChannelFit(:,1)=allChannelFit(:,1)*m(1).rangeScaleFac;
  for i=1:2
    %channelPrefs = e.stimVals;
    xp= channelPrefs(1):1:channelPrefs(end)*m(1).rangeScaleFac;
    yp=vonMises(allChannelFit(i,1:4),d2r(xp));
    xpOri=xp/m(1).rangeScaleFac;
    plot(channelPrefs(1:M), reshape(totalChannelTuning(1,1,i,:),[8,1]),[getcolor(i),'o']); hold on
    %plot(xpOri,yp,[getcolor(i),'--']);
  end
  title('Channel tuning functions');

end

figure;
subplot(1,3,1);
plot(1:2,rho,'ro');
hold on;
plot(1:2,mean(rho,2),'bo');
title('rho');
axis([0,3,0,0.5]);

subplot(1,3,2);
plot(1:2,sigma,'ro');
hold on;
plot(1:2,mean(sigma,2),'bo');
title('sigma');
axis([0,3,0,1]);

subplot(1,3,3);
plot(1:2,tao,'ro');
hold on;
plot(1:2,mean(tao,2),'bo');
title('tao');
axis([0,3,0,2]);

mean_posterior = zeros(nStimList,nROI,2,M);
std_posterior = zeros(nStimList,nROI,2,M);
for stimNum = 1:nStimList
	for ROI = 1:nROI
		for condi = 1:2
			for stimShown = 1:M
				startIndex = 1 + (stimShown-1)*trialPerStim;
				endIndex = (stimShown)*trialPerStim;
				mean_posterior(stimNum,ROI,condi,stimShown) = circ_mean(reshape(posterior_mean(stimNum,ROI,condi,:,startIndex:endIndex),[nsubj*trialPerStim,1]))/180*(2*pi)/(2*pi)*180;
				std_posterior(stimNum,ROI,condi,stimShown) = mean(reshape(posterior_std(stimNum,ROI,condi,:,startIndex:endIndex),[nsubj*trialPerStim,1]));
			end
		end
	end
end


thisROI = 2;
lowHighInd = 1;
medianList = zeros(8,180);
plotIndex = 1;
figure;
for i = 1:8
	startIndex = 1 + (i-1)*trialPerStim;
	endIndex = (i)*trialPerStim;
	trialPosterior = reshape(posterior(1,thisROI,lowHighInd,:,startIndex:endIndex,:),[nsubj*trialPerStim,180]);
	medianList(i,:) = median(trialPosterior,1);
	if(i == plotIndex)

		%plot(trialPosterior',':');
		hold on;
		plot(medianList(i,:))
	end
end

%rotate through to 90;
roundedStimVals = round(e.stimVals-0.1)+1;
for i = 1:M
	startIndex = 1 + (i-1)*trialPerStim;
	endIndex = (i)*trialPerStim;
	if(i == 1)
		trialPosterior = circshift(reshape(posterior(1,thisROI,lowHighInd,:,startIndex:endIndex,:),[nsubj*trialPerStim,180]),90 - roundedStimVals(i),2);
	else
		trialPosterior = [trialPosterior; reshape(posterior(1,thisROI,lowHighInd,:,startIndex:endIndex,:),[nsubj*trialPerStim,180])];
	end
end
medianList(i,:) = median(trialPosterior,1);
if(plotIndex == i)
	figure;
	plot(trialPosterior,':');
	hold on;
	plot(medianList(i,:))
end
