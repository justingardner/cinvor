% testCinvor.m
%
%      usage: testCinvor()
%         by: justin gardner
%       date: 09/07/16
%    purpose: function to test build/test channels based on cinvor data
%
function retval = testCinvor(varargin)

% get arguments
getArgs(varargin,{'dataDir=~/data/cinvor2','analName=decon1gIns','nFold=5','subjectList',{'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'},'recompute=0','figDir=~/Desktop'});

% number of subjects
nSubjects = length(subjectList);


dispFig3(dataDir);


% runf data analysis
contrastNames = {'high','low'};
for iContrast = 1:length(contrastNames)
  % get contrast name
  contrastName = contrastNames{iContrast};

  % get data
  [leftChannels rightChannels leftOutputs rightOutputs roiName nROI] = getCinvorData(dataDir,contrastName,analName,nFold,subjectList,nSubjects,recompute);

  % do subject-by-subject fit
  subjectBySubjectFit(iContrast) = cinvorSubjectBySubject(leftChannels,rightChannels,leftOutputs,rightOutputs);
  
  % display data
  [contraCombinedOutput(iContrast) contraCombinedChannel(iContrast) ipsiCombinedOutput(iContrast) ipsiCombinedChannel(iContrast)] = dispCinvorData(leftChannels,rightChannels,leftOutputs,rightOutputs,contrastName,roiName,nROI,nSubjects);

end

% now fit high-contrast and low-contrast data
[~,highContrastResponse,highContrastSTE,xVals] = dispChannelOutput(contraCombinedOutput(1),contraCombinedChannel(1),'suppressPlot=1');
[~,lowContrastResponse,lowContrastSTE,xVals] = dispChannelOutput(contraCombinedOutput(2),contraCombinedChannel(2),'suppressPlot=1');

% fit high / low contrast
highFit = fitVonMises(xVals,highContrastResponse,'dispFit=0');
lowFit = fitVonMises(xVals,lowContrastResponse,'dispFit=0');

% colors
colors.contraHigh = [1 0 0];
colors.contraLow = [242 149 9]/255;
colors.ipsiHigh = [0 0 0];
colors.ipsiLow = [0 0 0];

% display figure 1
dispFig1(figDir,contraCombinedOutput,contraCombinedChannel,ipsiCombinedOutput,ipsiCombinedChannel,highFit,lowFit,highContrastSTE,lowContrastSTE,colors);

% display figure 6
dispFig6(figDir,contraCombinedOutput,contraCombinedChannel,ipsiCombinedOutput,ipsiCombinedChannel,colors);

% statistics
disp(sprintf('Contra high contrast r2: %0.2f +- %0.2f',contraCombinedOutput(1).r2.overall,contraCombinedOutput(1).r2.steOverall));
disp(sprintf('Ipsi high contrast r2: %0.2f +- %0.2f',ipsiCombinedOutput(1).r2.overall,ipsiCombinedOutput(1).r2.steOverall));
disp(sprintf('Contra low contrast r2: %0.2f +- %0.2f',contraCombinedOutput(2).r2.overall,contraCombinedOutput(2).r2.steOverall));
disp(sprintf('Ipsi low contrast r2: %0.2f +- %0.2f',ipsiCombinedOutput(2).r2.overall,ipsiCombinedOutput(2).r2.steOverall));

%%%%%%%%%%%%%%%%%%
%    dispFig6    %
%%%%%%%%%%%%%%%%%%
function dispFig6(figDir,contraCombinedOutput,contraCombinedChannel,ipsiCombinedOutput,ipsiCombinedChannel,colors)

mlrSmartfig('cinvor_fig6','reuse');
clf;hold on

subplot(1,2,1);
dispChannelOutput(contraCombinedOutput(1),contraCombinedChannel(1),'likelihood=1','MarkerFaceColor',colors.ipsiHigh);
dispChannelOutput(ipsiCombinedOutput(1),ipsiCombinedChannel(1),'likelihood=1','MarkerFaceColor',colors.contraHigh);
title(sprintf('V1 high contrast'));
ylabel(sprintf('Decoded likelihood\np(orientation|BOLD)'));
xlabel('Orientation difference from true (deg)');
drawPublishAxis('xTick',[-90:90:90],'xMinorTick',[-90:45:90],'whichAxis=horizontal');

subplot(1,2,2);
dispChannelOutput(ipsiCombinedOutput(2),ipsiCombinedChannel(2),'likelihood=1','MarkerFaceColor',colors.ipsiLow);
dispChannelOutput(contraCombinedOutput(2),contraCombinedChannel(2),'likelihood=1','MarkerFaceColor',colors.contraLow);
title(sprintf('V1 low contrast'));
xlabel('Orientation difference from true (deg)');
drawPublishAxis('xTick',[-90:90:90],'xMinorTick',[-90:45:90],'whichAxis=horizontal','figSize=1');

% print figure
drawnow
print('-dpdf',fullfile(figDir,'fig6.pdf'));

%%%%%%%%%%%%%%%%%%
%    dispFig1    %
%%%%%%%%%%%%%%%%%%
function dispFig1(figDir,contraCombinedOutput,contraCombinedChannel,ipsiCombinedOutput,ipsiCombinedChannel,highFit,lowFit,highContrastSTE,lowContrastSTE,colors);


% plot
mlrSmartfig('cinvor_fig1','reuse');
clf;hold on
subplot(1,3,1);
dispChannelOutput(contraCombinedOutput(1),contraCombinedChannel(1),'MarkerFaceColor',colors.contraHigh);
dispChannelOutput(ipsiCombinedOutput(1),ipsiCombinedChannel(1),'MarkerFaceColor',colors.ipsiHigh);
title(sprintf('V1 high contrast'));
ylabel('Channel Response (percentile of full)');
drawPublishAxis('xTick',[-90:90:90],'xMinorTick',[-90:45:90],'yTick',[0:0.1:0.4]);

subplot(1,3,2);
dispChannelOutput(contraCombinedOutput(2),contraCombinedChannel(2),'MarkerFaceColor',colors.contraLow);
dispChannelOutput(ipsiCombinedOutput(2),ipsiCombinedChannel(2),'MarkerFaceColor',colors.ipsiLow);
title(sprintf('V1 low contrast'));
xlabel('Orientation difference from true (deg)');
drawPublishAxis('xTick',[-90:90:90],'xMinorTick',[-90:45:90],'yTick',[0:0.1:0.4],'whichAxis=horizontal');

subplot(1,3,3);
myerrorbar(lowFit.x,lowFit.y,'yError',lowContrastSTE,'Symbol','ko','MarkerFaceColor',colors.contraLow,'MarkerEdgeColor','w','MarkerSize',8.5)
plot(lowFit.xSmooth,lowFit.yFitSmooth,'k-','LineWidth',2,'Color',colors.contraLow);
myerrorbar(highFit.x,highFit.y,'yError',highContrastSTE,'Symbol','ro','MarkerFaceColor',colors.contraHigh,'MarkerEdgeColor','w','MarkerSize',8.5)
plot(highFit.xSmooth,highFit.yFitSmooth,'r-','LineWidth',2,'Color',colors.contraHigh);
drawPublishAxis('xTick',[-90:90:90],'xMinorTick',[-90:45:90],'yTick',[0:0.1:0.4],'whichAxis=horizontal','figSize=2');

mylegend({'Contra high contrast','Contra low contrast','Ipsi'},{{'ro' 'MarkerFaceColor' colors.contraHigh 'MarkerEdgeColor=w'},{'ro' 'MarkerFaceColor' colors.contraLow 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});

% print figure
drawnow
print('-dpdf',fullfile(figDir,'fig1.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispChannelLikelihood    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispChannelLikelihood(leftOutput,leftChannel,rightOutput,rightChannel,contrastName,subjectName,roiName,titleStr,figStr,nCol,iCol,contraIpsi)

if nargin < 10
  nCol = 1;
  iCol = 1;
end
if nargin < 12
  contraIpsi = false;
end

mlrSmartfig(figStr,'reuse');
if iCol == 1,clf;end

subplot(2,nCol,iCol);
dispChannelOutput(leftOutput,leftChannel);
dispChannelOutput(rightOutput,rightChannel,'MarkerFaceColor=r');
title(sprintf('%s: %s (%s)',roiName,subjectName,titleStr));
xlabel('Orientation difference from true (deg)');ylabel('Channel Response (percentile of full)');
if ~contraIpsi
  mylegend({'Right visual field','Left visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
else
  mylegend({'Ipsi visual field','Contra visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
end
drawPublishAxis('xTick',[-90:90:90],'xMinorTick',[-90:45:90],'yTick',[0:0.1:0.4]);

subplot(2,nCol,iCol+nCol);
statsStrLeft = dispChannelOutput(leftOutput,leftChannel,'likelihood=1');
statsStrRight = dispChannelOutput(rightOutput,rightChannel,'MarkerFaceColor=r','likelihood=1');
title(sprintf('Likelihood\n (left: %s)\n(right: %s)',statsStrLeft,statsStrRight));
xlabel('Orientation difference from true (deg)');ylabel('Likelihood (p)');

if ~contraIpsi
  mylegend({'Right visual field','Left visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
else
  mylegend({'Ipsi visual field','Contra visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
end

drawPublishAxis('xTick',[-90:90:90],'xMinorTick',[-90:45:90]);
drawnow

%%%%%%%%%%%%%%%%%%%%%%%
%    getCinvorData    %
%%%%%%%%%%%%%%%%%%%%%%%
function [leftChannels rightChannels leftOutputs rightOutputs roiName nROI] = getCinvorData(dataDir,contrastName,analName,nFold,subjectList,nSubjects,recompute)

% save name
saveName = fullfile(dataDir,sprintf('testCinvorFull%s',contrastName));

% if we can reload then reload
if ~recompute && isfile(setext(saveName,'mat'))
  % just load previously computed
  disp(sprintf('(testCinvor) Loading precomputed analysis: %s',saveName));
  load(saveName);
else
  % display what we are doing
  disppercent(-inf,sprintf('(testCinvor) Building and testing channels for %s contrast using analysis %s with %i fold cross-validation on %i subjects',contrastName,analName,nFold,nSubjects));

  % iterate over subjects
  for iSubject = 1:nSubjects
    % load the data
    dataName = fullfile(dataDir,subjectList{iSubject},'Anal',analName);
    dispHeader
    disp(sprintf('(testCinvor) Loading data for subject %i/%i: %s',iSubject,nSubjects,dataName));
    dispHeader
    load(dataName);
  
    % specify the experiment
    e = specifyCinvorExperiment('stimLevel=8','trialPerStim=1');
 
    % get number of ROIS
    nROI = length(s.lvf.roi);

    % decide on which contrast to do
    if strcmp(lower(contrastName),'high');
      stimNums = e.highContrastStimuli;
    else
      stimNums = e.lowContrastStimuli;
    end
    
    % cycle over ROIs
    for iROI = 1:nROI
      
      % grab instance for this ROI left and right visual field
      roiName{iROI} = s.lvf.roi{iROI}.name;
      stimNames{iROI} = {s.lvf.stimNames{stimNums}};
      lvfInstances = {s.lvf.roi{iROI}.instance.instances{stimNums}};
      rvfInstances = {s.rvf.roi{iROI}.instance.instances{stimNums}};
      
      % get the crossVal sets
      crossVal = getCrossValSets(lvfInstances,'nFold',nFold);
      
      % now train and test over folds
      for iFold = 1:nFold
	% comupte train/test of left visual field channels
	[trainInstances testInstances] = getCrossValInstances(lvfInstances,crossVal,iFold);
	leftFoldChannel(iFold) = buildChannels(trainInstances,e.stimVals,'fitNoiseModel=1','noiseModelGridSearchOnly=0','noiseModelFitTolerence',0.1,'noiseModelGridSteps=25');
	leftFoldOutputs(iFold) = testChannels(testInstances,e.stimVals,leftFoldChannel(iFold));
	

	% comupte train/test of right visual field channels
	[trainInstances testInstances] = getCrossValInstances(rvfInstances,crossVal,iFold);
	rightFoldChannel(iFold) = buildChannels(trainInstances,e.stimVals,'fitNoiseModel=1','noiseModelGridSearchOnly=0','noiseModelFitTolerence',0.1,'noiseModelGridSteps=25');
	rightFoldOutputs(iFold) = testChannels(testInstances,e.stimVals,rightFoldChannel(iFold));

	% display
	dispChannelLikelihood(leftFoldOutputs(iFold),leftFoldChannel(iFold),rightFoldOutputs(iFold),rightFoldChannel(iFold),contrastName,subjectList{iSubject},roiName{iROI},sprintf('Fold: %i/%i',iFold,nFold),sprintf('%s%s%sFold',contrastName,subjectList{iSubject},roiName{iROI}),nFold,iFold);

	% disp percent
	disppercent(calcPercentDone(iSubject,nSubjects,iROI,nROI,iFold,nFold));
      end

      % average over folds and store in structure
      [leftOutputs(iSubject,iROI) leftChannels(iSubject,iROI)] = combineChannelOutputs(leftFoldOutputs,leftFoldChannel);
      [rightOutputs(iSubject,iROI) rightChannels(iSubject,iROI)] = combineChannelOutputs(rightFoldOutputs,rightFoldChannel);
      % display
      dispChannelLikelihood(leftOutputs(iSubject,iROI),leftChannels(iSubject,iROI),rightOutputs(iSubject,iROI),rightChannels(iSubject,iROI),contrastName,subjectList{iSubject},roiName{iROI},sprintf('Across folds'),sprintf('%s%s',contrastName,subjectList{iSubject}),nROI,iROI);
    end
  end

  disppercent(inf);
  disp(sprintf('(testCinvor) Saving precompute file: %s',saveName));
  save(saveName);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    dispCinvorData    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [contraCombined contraCombinedChannel ipsiCombined ipsiCombinedChannel] = dispCinvorData(leftChannels,rightChannels,leftOutputs,rightOutputs,contrastName,roiName,nROI,nSubjects)

for iROI = 1:2:nROI
  % get contras/ipsi for left
  [contraLeft contraLeftChannel] = combineChannelOutputs(rightOutputs(:,iROI),rightChannels(:,iROI));
  [ipsiLeft ipsiLeftChannel] = combineChannelOutputs(leftOutputs(:,iROI),leftChannels(:,iROI));

  % get contras/ipsi for right
  [contraRight contraRightChannel] = combineChannelOutputs(leftOutputs(:,iROI+1),leftChannels(:,iROI+1));
  [ipsiRight ipsiRightChannel]= combineChannelOutputs(rightOutputs(:,iROI+1),rightChannels(:,iROI+1));

  % get combined contra / ipsi
  [contraCombined contraCombinedChannel] = combineChannelOutputs([rightOutputs(:,iROI) leftOutputs(:,iROI+1)],[rightChannels(:,iROI) leftChannels(:,iROI+1)]);
  [ipsiCombined ipsiCombinedChannel] = combineChannelOutputs([rightOutputs(:,iROI+1) leftOutputs(:,iROI)],[leftChannels(:,iROI+1) rightChannels(:,iROI)]);;

  % display
  dispChannelLikelihood(contraLeft,contraLeftChannel,ipsiLeft,ipsiLeftChannel,contrastName,sprintf('Average (n=%i)',nSubjects),roiName{iROI},'left',sprintf('%s_average',contrastName),nROI+nROI/2,iROI);
  dispChannelLikelihood(ipsiRight,ipsiRightChannel,contraRight,contraRightChannel,contrastName,sprintf('Average (n=%i)',nSubjects),roiName{iROI},'right',sprintf('%s_average',contrastName),nROI+nROI/2,iROI+1);
  dispChannelLikelihood(contraCombined,contraCombinedChannel,ipsiCombined,ipsiCombinedChannel,contrastName,sprintf('Average (n=%i)',nSubjects),roiName{iROI},'Combined left+right',sprintf('%s_average',contrastName),nROI+nROI/2,iROI+2,true);

  % make equal axis
  makeEqualYaxis(2,nROI+nROI/2,iROI:iROI+2);
  makeEqualYaxis(2,nROI+nROI/2,nROI+nROI/2+(iROI:iROI+2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    cinvorSubjectBySubject    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = cinvorSubjectBySubject(leftChannels,rightChannels,leftOutputs,rightOutputs)

% go subject-by-subjct
for iSubject = 1:size(leftChannels,1)

  % collect left contralateral
  [~,contraResponse(2*iSubject-1,:),contrastResponseSTE(2*iSubject-1,:),xVals] = dispChannelOutput(leftOutputs(iSubject,2),leftChannels(iSubject,2),'suppressPlot=1');
  
  % collect right contralateral
  [~,contraResponse(2*iSubject,:),contrastResponseSTE(2*iSubject,:),xVals] = dispChannelOutput(rightOutputs(iSubject,1),rightChannels(iSubject,1),'suppressPlot=1');
  
  % collect left ipsilateral
  [~,ipsiResponse(2*iSubject-1,:),ipsiResponseSTE(2*iSubject-1,:),xVals] = dispChannelOutput(leftOutputs(iSubject,1),leftChannels(iSubject,1),'suppressPlot=1');
  
  % collect right ipsilateral
  [~,ipsiResponse(2*iSubject,:),ipsiResponseSTE(2*iSubject,:),xVals] = dispChannelOutput(rightOutputs(iSubject,2),rightChannels(iSubject,2),'suppressPlot=1');
  
end

disppercent(-inf,'(testCinvor:cinvorSubjectBySubject) Computing von mises fits for each hemisphere');
for iHemisphere = 1:size(contraResponse,1)
  retval.fit(iHemisphere) = fitVonMises(xVals,contraResponse(iHemisphere,:),'dispfit=0');
  retval.kappa(iHemisphere) = retval.fit(iHemisphere).params.kappa;
  retval.mu(iHemisphere) = retval.fit(iHemisphere).params.mu;
  retval.offset(iHemisphere) = retval.fit(iHemisphere).params.offset;
  retval.amp(iHemisphere) = retval.fit(iHemisphere).params.amp;
  retval.kappaVar(iHemisphere) = retval.fit(iHemisphere).covar(2,2);
  disppercent(iHemisphere/size(contraResponse,1));
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    simulateNeuralTuningWidths    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simResults = simulateNeuralTuningWidths(dataDir,varargin)

% get default arguments
filterType = [];channelExponent = [];nNoiseVals = [];n = [];filterType = [];recompute = [];
getArgs(varargin,{'n=1000','nNoiseVals=50','channelExponent=7','filterType=sinFilter','recompute=0','nStimVals=8','numFilters=[]'});

% compute number of filters needed for exponent
if isempty(numFilters)
  numFilters = channelExponent+1;
end

% set the experiment (number of orientations to test and repetitions)
simResults.e = specifyCinvorExperiment('stimLevel',nStimVals,'trialPerStim=21');
simResults.numFilters = numFilters;
simResults.channelExponent = channelExponent;
simResults.filterType = filterType;

% Setting kappa values to use
% using fitVonMises to convert from half-width-at-half-height in degrees to kappa
simResults.kappa = fitVonMises([],[],'halfWidthAtHalfHeight',15:5:60);
nKappa = length(simResults.kappa);

% values of noise on voxels
simResults.noiseVals = [0 logspace(log10(0.025),log10(15),nNoiseVals-1)];

% make simulation name to save under
simName = sprintf('kappaSimulation_n%i_noise_%i_e%i_%s.mat',n,nNoiseVals,channelExponent,filterType);
simName = fullfile(dataDir,simName);


if ~isequal(recompute,1) && isfile(simName)
  disp(sprintf('(testCinvor) Loading precomputed file: %s',simName));
  load(simName);
else
  % get / set number of workers for the parfor below
  disp(sprintf('(testCinvor) Number of parallel workers: %i',mlrNumWorkers(1)));
  % recompute
  dispHeader(sprintf('(testCinvor:simulateNeuralTuningWdiths)\nExponent: %i filterType: %s nStimVals: %i\nNoise levels: %i Kappa values: %i Num simulations: %i',channelExponent,filterType,nStimVals,nNoiseVals,nKappa,n));
  disppercent(-inf,sprintf('%i total simulations',nNoiseVals*nKappa*n));

  % init variables
  simResults.r2 = nan(nKappa,nNoiseVals,n);
  halfWidth = nan(1,n);
  simResults.halfWidth = nan(nKappa,nNoiseVals,n);
  % over noise values
  for iNoise = 1:nNoiseVals
    % run simulation over tuning widths
    for iKappa = 1:nKappa
      % set the model parameters
      m = setCinvorModel(simResults.e,'kappa',simResults.kappa(iKappa),'noise',simResults.noiseVals(iNoise));
      % n simulations
      parfor iSimulation = 1:n
	warning off;
	% get train and test instances
	trainInstances = getCinvorInstances(m,simResults.e);
	testInstances = getCinvorInstances(m,simResults.e);
	% compute forward model
	channel(iSimulation) = buildChannels(trainInstances,simResults.e.stimVals,'dispChannels=0','fitNoise',0,'model',filterType,'exponent',channelExponent,'numFilters',numFilters);
	% test on test instances
	channelOutput(iSimulation) = testChannels(testInstances,simResults.e.stimVals,channel(iSimulation),'fitNoise',0);
	% do von mises fit
	vonMisesFit(iSimulation) = fitVonMises(channelOutput(iSimulation).averageChannelResponseXvals,channelOutput(iSimulation).averageChannelResponse,'dispFit=0');
	halfWidth(iSimulation) = vonMisesFit(iSimulation).params.halfWidthAtHalfHeight;
      end
      % init channel responses if necessary
      nXvals = length(channelOutput(1).averageChannelResponseXvals);
      if ~isfield(simResults,'channelResponse')
	simResults.channelResponse = nan(nKappa,nNoiseVals,n,nXvals);
      end
      % save channel Responses
      simResults.channelResponse(iKappa,iNoise,:,:) = reshape([channelOutput(:).averageChannelResponse],nXvals,n)';
      simResults.vonMisesFit(iKappa,iNoise,:) = vonMisesFit;
      simResults.halfWidth(iKappa,iNoise,:) = halfWidth;
      simResults.r2(iKappa,iNoise,:) = [channelOutput(:).r2];
      % put other variables into siResults
      simResults.m(iKappa,iNoise) = m;
      simResults.channel = channel(1);
      simResults.channelXvals = channelOutput(1).averageChannelResponseXvals;
      % display done
      disppercent(calcPercentDone(iNoise,nNoiseVals,iKappa,nKappa));
      %channel,m,channelResponse,fit,r2,vonMisesFit,halfWidth
    end
  end
  disppercent(inf);
  save(simName,'simResults');
end
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispSimulateNeuralTuningWidths    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispSimulateNeuralTuningWidths(simResults,varargin)

getArgs(varargin,{'f1=[]','f2=[]','subplotInfo=[]'});

if isempty(f1)
  f1 = mlrSmartfig('testCinvorTuningWidths','reuse');
  clf(f1);
end
if isempty(f2)
  f2 = mlrSmartfig('testCinvorTuningWidths2','reuse');
  clf(f2);
end

% calculate actual tuning width
[~,stimIndex] = min(abs(simResults.channel.idealStimVals-90));
channelResponse = simResults.channel.idealStimResponse(stimIndex,:);
x = simResults.channel.channelPref;
fit = fitVonMises(x-90,channelResponse,'dispFit=0');
disp(sprintf('(testCinvor:dispSimulateNeuralTuningWidths) Actual tuning width: %0.2f',fit.params.halfWidthAtHalfHeight));

% get mean channel response
xVals = simResults.channelXvals;
meanResponse = squeeze(mean(simResults.channelResponse,3));

figure(f1);
for iKappa = 1:length(simResults.kappa)
  % grab r2 and halfWidth values
  r2 = squeeze(simResults.r2(iKappa,:,:));
  halfWidth = squeeze(simResults.halfWidth(iKappa,:,:));
  % bin by r2
  r2binWidth = 0.05;
  r2bins = [0:r2binWidth:(1-r2binWidth/2) 1];
  for iBin = 1:length(r2bins)
    % get the median value across the bin
    halfWidthBin(iBin) = median(halfWidth((r2(:)>=r2bins(iBin)) & (r2(:)<(r2bins(iBin)+r2binWidth))));
  end
  
  % color to plot width
  c = getSmoothColor(iKappa,length(simResults.kappa),'cool');
  % plot summaries
  if ~isempty(subplotInfo)
    subplot(subplotInfo(1),subplotInfo(2),subplotInfo(3));
  end
  % plot
  plot([r2bins(1:end-1)+r2binWidth/2 1],halfWidthBin,'o-','MarkerFaceColor',c,'MarkerEdgeColor','w','Color',c);hold on
end


% get tuning for different r2 values
targetr2 = [0.25 0.5 1];r2binwidth = 0.25;
figure(f2)
for iR2 = 1:length(targetr2)
  for iKappa = 1:length(simResults.kappa)
    % get all simulations with the desired r2 in the bin
    r2 = squeeze(simResults.r2(iKappa,:,:));
    [noiseVals,simulationNums] = find((r2 > (targetr2(iR2)-r2binwidth)) & (r2 < (targetr2(iR2)+r2binwidth)));
    %disp(sprintf('(testCinvor:dispSimulateNeuralTuningWidths) Found %i simulations for r2: %0.3f',length(noiseVals),targetr2(iR2)));
    % no make a matrix of all of their channel responses
    channelResponses = [];
    for iSim = 1:length(noiseVals)
      channelResponses(iSim,:) = squeeze(simResults.channelResponse(iKappa,noiseVals(iSim),simulationNums(iSim),:));
    end

    % color to plot width
    c = getSmoothColor(iKappa,length(simResults.kappa),'cool');
    % plot average
    % and plot average
    subplot(1,length(targetr2),iR2);
    plot(simResults.channelXvals,mean(channelResponses),'o-','MarkerEdgeColor','w','Color',c);hold on
  end
end

% label axis
figure(f1)
if ~isempty(subplotInfo)
  subplot(subplotInfo(1),subplotInfo(2),subplotInfo(3));
end
hline(fit.params.halfWidthAtHalfHeight);
xlabel('r2');  ylabel('Channel model width');
a = axis;
yaxis(0,40);
drawPublishAxis;

figure(f2)
for iSubplot = 1:length(targetr2)
  subplot(1,length(targetr2),iSubplot);
  xlabel('Orientation difference from true (deg)');
  ylabel('Channel response');
  title(sprintf('r2=%0.2f',targetr2(iSubplot)));
  drawPublishAxis('figSize=2')
end





%%%%%%%%%%%%%%%%%
%    dispFig3   %
%%%%%%%%%%%%%%%%%
function dispFig3(dataDir)

% simulate range of neural tuning widths
simResults3 = simulateNeuralTuningWidths(dataDir,'n=500','nNoiseVals=50','recompute=0','channelExponent=3','nStimVals=8','numFilters=8');
simResultsStick = simulateNeuralTuningWidths(dataDir,'n=500','nNoiseVals=50','recompute=0','channelExponent=7','nStimVals=8','numFilters=8','filterType=stickFilter');
simResults7 = simulateNeuralTuningWidths(dataDir,'n=500','nNoiseVals=50','recompute=0','channelExponent=7','nStimVals=8','numFilters=8');

%dispSimulateNeuralTuningWidths(simResults);
keyboard


f = mlrSmartfig('testCinvor3','reuse');clf(f);
dispSimulateNeuralTuningWidths(simResults3,'f1',f,'subplotInfo',[1 3 1]);
dispSimulateNeuralTuningWidths(simResults7,'f1',f,'subplotInfo',[1 3 2]);
dispSimulateNeuralTuningWidths(simResultsStick,'f1',f,'subplotInfo',[1 3 3]);
