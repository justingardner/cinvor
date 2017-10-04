% testCinvor.m
%
%      usage: testCinvor()
%         by: justin gardner
%       date: 09/07/16
%    purpose: function to test build/test channels based on cinvor data
%
function retval = testCinvor(varargin)

% get arguments
getArgs(varargin,{'dataDir=~/Google Drive/cinvor/precompute','analName=decon1gIns','nFold=5','subjectList',{'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'},'recompute=0','figDir=~/Desktop'});

% number of subjects
nSubjects = length(subjectList);

n = 100; nNoiseVals = 40; recompute = 0;
recompute = 0;n = 50;nNoiseVals=10;
simResults = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=7','nStimVals=8','numFilters=8','kappa=[-8]','filterType=stickFilter','minNoise=0.001','maxNoise=1','categoryConfusion=0.15');
dispSimulateNeuralTuningWidths(simResults);
keyboard
%dispSimulateNeuralTuningWidths(simResults,'r2plot=0');;
dispChannelTuning(simResults,0.1,0.1)
keyboard

%keyboard
% run simulations
%computeSimulations(dataDir,'n=100','nNoiseVals=40','recompute=0');
%simResults = computeSimulations(dataDir,'n=500','nNoiseVals=50','recompute=0');
simResults = computeSimulations(dataDir,'n=1000','nNoiseVals=50','recompute=0');

dispChannelWidthByR2Fig({simResults{1:3}},dataDir,figDir);
%dispR2ByTuningWidthFig(simResults{2},dataDir,figDir);

dispSimulateNeuralTuningWidths(simResults{7});

dispCategoryFigure({simResults{[5 4 6]}},dataDir,figDir);
dispChannelWidthByR2Fig({simResults{1:3}},dataDir,figDir);

% run data analysis
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

  % get r2 for subject-by-subject
  for iSubject = 1:size(leftOutputs,1)
    r2(iContrast,iSubject*2-1) = leftOutputs(iSubject,2).r2.overall;
    r2(iContrast,iSubject*2) = rightOutputs(iSubject,1).r2.overall;
  end
end

% now fit high-contrast and low-contrast data
[~,highContrastResponse,highContrastSTE,xVals] = dispChannelOutput(contraCombinedOutput(1),contraCombinedChannel(1),'suppressPlot=1');
[~,lowContrastResponse,lowContrastSTE,xVals] = dispChannelOutput(contraCombinedOutput(2),contraCombinedChannel(2),'suppressPlot=1');

% fit high / low contrast
highFit = fitVonMises(xVals,highContrastResponse,'dispFit=0','parametricBootstrap=1000','yste',highContrastSTE);
lowFit = fitVonMises(xVals,lowContrastResponse,'dispFit=0','parametricBootstrap=1000','yste',lowContrastSTE);

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

disp(sprintf('high contrast half-width: %0.1f',highFit.params.halfWidthAtHalfHeight));
disp(sprintf('low contrast half-width: %0.1f',lowFit.params.halfWidthAtHalfHeight));
% calculate percentage of bootstraps that have high narrower than low
p = sum((highFit.yBootHalfWidth - lowFit.yBootHalfWidth)<0)/length(lowFit.yBootHalfWidth);
disp(sprintf('Percent bootstraps with high narrower than low: %0.1f',p*100));


disp(sprintf('high contrast amplitude: %0.2f',highFit.params.amp));
disp(sprintf('low contrast amplitude: %0.2f',lowFit.params.amp));
p = sum((highFit.yBoot(:,1) - lowFit.yBoot(:,1))>0)/length(lowFit.yBootHalfWidth);
disp(sprintf('Percent bootstraps with high higher than low: %0.1f',p*100));

% compute hemisphere stats
lowContrastWider = 0;lowContrastLower = 0;
nHemi = length(subjectBySubjectFit(1).fit);
for iHemi = 1:nHemi
  % half width
  halfWidth(1,iHemi) = subjectBySubjectFit(1).fit(iHemi).params.halfWidthAtHalfHeight;
  halfWidth(2,iHemi) = subjectBySubjectFit(2).fit(iHemi).params.halfWidthAtHalfHeight;
  % amplitude
  amp(1,iHemi) = subjectBySubjectFit(1).fit(iHemi).params.amp;
  amp(2,iHemi) = subjectBySubjectFit(2).fit(iHemi).params.amp;
end
disp(sprintf('%i/%i hemispheres have wider fit for low contrast',sum(halfWidth(1,:) < halfWidth(2,:)),size(halfWidth,2)));
[p t dof] = pairedttest(halfWidth(1,:),halfWidth(2,:));
disp(sprintf('p = %0.3f, t(%i) = %0.3f paired t-test',p,dof,t));

disp(sprintf('%i/%i hemispheres have wider fit for low contrast',sum(amp(1,:) > amp(2,:)),size(amp,2)));
[p t dof] = pairedttest(amp(1,:),amp(2,:));
disp(sprintf('p = %0.3f, t(%i) = %0.3f paired t-test',p,dof,t))

disp(sprintf('%i/%i hemispheres have higher r2 for high contrast',sum(r2(1,:) > r2(2,:)),size(r2,2)));
[p t dof] = pairedttest(r2(1,:),r2(2,:));
disp(sprintf('p = %0.3f, t(%i) = %0.3f paired t-test',p,dof,t))

% number taken from Dylan's analysis - simuCinvorChangeWidth.m
disp(sprintf('Magnitude decrease from high contrast is %0.3f%%',100*0.73/1.73));


%%%%%%%%%%%%%%%%%%
%    dispFig6    %
%%%%%%%%%%%%%%%%%%
function dispFig6(figDir,contraCombinedOutput,contraCombinedChannel,ipsiCombinedOutput,ipsiCombinedChannel,colors)

mlrSmartfig('cinvor_fig6','reuse');
clf;hold on

subplot(1,2,1);
dispChannelOutput(contraCombinedOutput(1),contraCombinedChannel(1),'likelihood=1','MarkerFaceColor',colors.contraHigh);
dispChannelOutput(ipsiCombinedOutput(1),ipsiCombinedChannel(1),'likelihood=1','MarkerFaceColor',colors.ipsiHigh);
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
  retval.fit(iHemisphere) = fitVonMises(xVals,contraResponse(iHemisphere,:),'dispFit=0');
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
getArgs(varargin,{'n=1000','nNoiseVals=50','channelExponent=7','filterType=sinFilter','recompute=0','nStimVals=8','numFilters=[]','weighting=random','kConcentrated=[]','kappa=[]','nNeuronTypes=180','minNoise=0.001','maxNoise=1','categoryConfusion=0'});

% compute number of filters needed for exponent
if isempty(numFilters)
  numFilters = channelExponent+1;
end

% set the experiment (number of orientations to test and repetitions)
simResults.e = specifyCinvorExperiment('stimLevel',nStimVals,'trialPerStim=21');
simResults.numFilters = numFilters;
simResults.channelExponent = channelExponent;
simResults.filterType = filterType;
simResults.kConcentrated = kConcentrated;

% Setting kappa values to use
% using fitVonMises to convert from half-width-at-half-height in degrees to kappa
if isempty(kappa)
  % set range of kappas
  simResults.kappa = [inf fitVonMises([],[],'halfWidthAtHalfHeight',5:5:45)];
  simName = 'kappa';
else
  % use passed in kappas
  simResults.kappa = kappa;
  simName = sprintf('kappa_%s',fixBadChars(num2str(kappa)));
end
nKappa = length(simResults.kappa);

% values of noise on voxels
simResults.noiseVals = [0 logspace(log10(minNoise),log10(maxNoise),nNoiseVals-1)];

% make simulation name to save under
dispStr = 'Random uniform weighting';
simName = sprintf('%sSimulation_n%i_noise_%i_e%i_%s',simName,n,nNoiseVals,channelExponent,filterType);
if ~isempty(kConcentrated) && strcmp(lower(weighting),'concentrated')
  simName = sprintf('%s_kConcentrate_%ideg',simName,round(fitVonMises([],[],'kappa',kConcentrated)));
  dispStr = sprintf('Concentrated weighting: %i',round(fitVonMises([],[],'kappa',kConcentrated)));
end
if ~isequal(categoryConfusion,0)
  simName = sprintf('%s_catcon_%i%',simName,round(100*categoryConfusion));
end
simName = fullfile(dataDir,simName);
simName = setext(simName,'mat');

if ~isequal(recompute,1) && isfile(simName)
  disp(sprintf('(testCinvor) Loading precomputed file: %s',simName));
  load(simName);
else
  % get / set number of workers for the parfor below
  disp(sprintf('(testCinvor) Number of parallel workers: %i',mlrNumWorkers(1)));
  % recompute
  dispHeader(sprintf('(testCinvor:simulateNeuralTuningWdiths)\nExponent: %i filterType: %s nStimVals: %i\nNoise levels: %i Kappa values: %i Num simulations: %i\n%s',channelExponent,filterType,nStimVals,nNoiseVals,nKappa,n,dispStr));

  disppercent(-inf,sprintf('%i total simulations',nNoiseVals*nKappa*n));

  % init variables
  simResults.r2 = nan(nKappa,nNoiseVals,n);
  halfWidth = nan(1,n);
  simResults.halfWidth = nan(nKappa,nNoiseVals,n);
  % over noise values
  for iNoise = 1:nNoiseVals
    % run simulation over tuning widths
    for iKappa = 1:nKappa
      % if kappa is less than zero, then we are doing n category tuning - so
      % set nTuning and neuronsPervox to the number of categories
      if simResults.kappa(iKappa) < 0
	nNeuronTypes = abs(simResults.kappa(iKappa));
      else
	% otherwise one neuron type for each orientation
	nNeuronTypes = 180;
      end
      % set the model parameters
      m = setCinvorModel(simResults.e,'kappa',simResults.kappa(iKappa),'noise',simResults.noiseVals(iNoise),'weighting',weighting,'kConcentrated',kConcentrated,'nTuning',nNeuronTypes,'neuronsPerVox',nNeuronTypes,'categoryConfusion',categoryConfusion);
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
      %simResults.vonMisesFit(iKappa,iNoise,:) = vonMisesFit;
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

getArgs(varargin,{'f1=[]','f2=[]','subplotInfo=[]','r2plot',1});

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
  if r2plot
    % grab r2 and halfWidth values
    r2 = squeeze(simResults.r2(iKappa,:,:));
    halfWidth = squeeze(simResults.halfWidth(iKappa,:,:));
    % bin by r2
    binWidth = 0.1;
    bins = [0:binWidth:(1-binWidth/2) 1];
    for iBin = 1:length(bins)
      % get the median value across the bin
      halfWidthBin(iBin) = median(halfWidth((r2(:)>=bins(iBin)) & (r2(:)<(bins(iBin)+binWidth))));
    end
    plotBins = bins;
  else
    bins = simResults.noiseVals;
    binWidth = 0;
    for iNoiseBin = 1:length(bins)
      % get the median value across the bin
      halfWidthBin(iNoiseBin) = median(simResults.halfWidth(iKappa,iNoiseBin,:));
    end
    plotBins = [bins(1:end-1)+binWidth/2 1];
  end
  % color to plot width
  c = getSmoothColor(iKappa,length(simResults.kappa),'cool');
  % plot summaries
  if ~isempty(subplotInfo)
    subplot(subplotInfo(1),subplotInfo(2),subplotInfo(3));
  end
  % plot
  plot(plotBins,halfWidthBin,'o-','MarkerFaceColor',c,'MarkerEdgeColor','w','Color',c);hold on
end

% get tuning for different r2 values
targetr2 = [0.25 0.5 1];r2binwidth = 0.25;
figure(f2);
% plot channel tuning graph for different r2
for iR2 = 1:length(targetr2)
  subplot(1,length(targetr2),iR2);
  dispChannelTuning(simResults,targetr2(iR2),r2binwidth)
end

% label axis
figure(f1)
if ~isempty(subplotInfo)
  subplot(subplotInfo(1),subplotInfo(2),subplotInfo(3));
end
hline(fit.params.halfWidthAtHalfHeight);
if r2plot
  xlabel('r2');
else
  xlabel('noise (std)');
  semilogx
end

ylabel('Channel model width');
a = axis;
yaxis(0,50);
drawPublishAxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispChannelTuning    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispChannelTuning(simResults,targetr2,r2binwidth)

for iKappa = 1:length(simResults.kappa)
  % get all simulations with the desired r2 in the bin
  r2 = squeeze(simResults.r2(iKappa,:,:));
  [noiseVals,simulationNums] = find((r2 > (targetr2-r2binwidth)) & (r2 < (targetr2+r2binwidth)));
  % no make a matrix of all of their channel responses
  channelResponses = [];
  for iSim = 1:length(noiseVals)
    channelResponses(iSim,:) = squeeze(simResults.channelResponse(iKappa,noiseVals(iSim),simulationNums(iSim),:));
  end

  % color to plot width
  c = getSmoothColor(iKappa,length(simResults.kappa),'cool');
  % plot average
  % and plot average
  plot(simResults.channelXvals,mean(channelResponses),'o-','MarkerEdgeColor','w','Color',c);hold on
end

xlabel('Orientation difference from true (deg)');
ylabel('Channel response');
title(sprintf('r2=%0.2f',targetr2));
drawPublishAxis('figSize=2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispR2ByTuningWidthFig    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispR2ByTuningWidthFig(simResults,dataDir,figDir)

mlrSmartfig('testCinvor_fig4','reuse');clf;

% get neural tuningwidth
neuralTuningWidth = fitVonMises([],[],'kappa',simResults.kappa);

for iNoise = 1:length(simResults.noiseVals)
  r2 = [];
  % get current color
  c = getSmoothColor(iNoise,length(simResults.noiseVals),'parula');
  for iKappa = 1:length(simResults.kappa)
    % get mean r1 value
    r2(iKappa) = median(simResults.r2(iKappa,iNoise,:));
  end
  plot(neuralTuningWidth,r2,'o-','MarkerFaceColor',c,'MarkerEdgeColor','w','Color',c);hold on
  legendNames{iNoise} = sprintf('%0.4f',simResults.noiseVals(iNoise));
  legendSymbols{iNoise} = {'o-' c 'MarkerFaceColor' c 'MarkerEdgeColor' 'w'};
end
xlabel('Neural tuning width (deg)');
ylabel('r2');
title('r2 for different SNR as a function of neural tuning width');
%mylegend(legendNames,legendSymbols);
h = colorbar;
set(h,'Limits',[min(simResults.noiseVals) max(simResults.noiseVals)]);
colorbarTicks = [0.001 0.01 0.1 1];
for iColorbarTicks = 1:length(colorbarTicks)
  colorbarTickLinear(iColorbarTicks) = find(colorbarTicks(iColorbarTicks)==simResults.noiseVals)/length(simResults.noiseVals);
  colorbarTickLabels{iColorbarTicks} = num2str(colorbarTicks(iColorbarTicks));
end
set(h,'Ticks',colorbarTickLinear)
set(h,'TickLabels',colorbarTickLabels)
yaxis(0,1);
drawPublishAxis;
print('-dpdf',fullfile(figDir,'fig4.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeSimulations    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simResults = computeSimulations(dataDir,varargin)

getArgs(varargin,{'n=500','nNoiseVals=50','recompute=0'});
		  
% simulate range of neural tuning widths
simResults{1} = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=3','nStimVals=8','numFilters=8');
simResults{1}.modelName = 'Channel exponent = 3';
simResults{2} = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=7','nStimVals=8','numFilters=8');
simResults{2}.modelName = 'Channel exponent = 7';
simResults{3} = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=7','nStimVals=8','numFilters=8','filterType=stickFilter');
simResults{3}.modelName = 'Stick filter';

% category model
simResults{4} = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=7','nStimVals=8','numFilters=8','kappa=[-2 -4 -8]');
simResults{4}.modelName = 'category, channel exponent = 7';
simResults{5} = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=3','nStimVals=8','numFilters=8','kappa=[-2 -4 -8]');
simResults{5}.modelName = 'category, stick filter';
simResults{6} = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=7','nStimVals=8','numFilters=8','kappa=[-2 -4 -8]','filterType=stickFilter');
simResults{6}.modelName = 'category';

% Compute models with weighting
simResults{7}= simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=7','nStimVals=8','numFilters=8','weighting=concentrated','kConcentrated',fitVonMises([],[],'halfWidthAtHalfHeight',30));
simResults{7}.modelName = 'Concentrated tuning 30';
simResults{8} = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=7','nStimVals=8','numFilters=8','weighting=concentrated','kConcentrated',fitVonMises([],[],'halfWidthAtHalfHeight',60));
simResults{8}.modelName = 'Concentrated tuning 60';
simResults{9} = simulateNeuralTuningWidths(dataDir,'n',n,'nNoiseVals',nNoiseVals,'recompute',recompute,'channelExponent=7','nStimVals=8','numFilters=8','weighting=concentrated','kConcentrated',inf);
simResults{9}.modelName = 'Concentrated tuning stick';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispChannelWidthByR2Fig   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispChannelWidthByR2Fig(simResults,dataDir,figDir)

f = mlrSmartfig('testCinvor3','reuse');clf(f);
subplot(2,3,1);
plot(simResults{1}.channel.channelPref-90,simResults{1}.channel.idealStimResponse(5,:),'ko-','MarkerFaceColor','k','MarkerEdgeColor','w');
xlabel('Channel pref (deg)');
ylabel('Channel response');
drawPublishAxis;

subplot(2,3,2);
plot(simResults{2}.channel.channelPref-90,simResults{2}.channel.idealStimResponse(5,:),'ko-','MarkerFaceColor','k','MarkerEdgeColor','w');
xlabel('Channel pref (deg)');
ylabel('Channel response');
drawPublishAxis;

subplot(2,3,3);
for iStick = 1:length(simResults{3}.channel.channelPref)
  plot([simResults{3}.channel.channelPref(iStick)-90 simResults{3}.channel.channelPref(iStick)-90],[0 simResults{3}.channel.idealStimResponse(5,iStick)],'k-','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',7);hold on
  plot(simResults{3}.channel.channelPref(iStick)-90,simResults{3}.channel.idealStimResponse(5,iStick),'ko','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',7);hold on
end
xlabel('Channel pref (deg)');
ylabel('Channel response');
drawPublishAxis;

dispSimulateNeuralTuningWidths(simResults{1},'f1',f,'subplotInfo',[2 3 4]);
dispSimulateNeuralTuningWidths(simResults{2},'f1',f,'subplotInfo',[2 3 5]);
dispSimulateNeuralTuningWidths(simResults{3},'f1',f,'subplotInfo',[2 3 6]);
neuralTuningWidth = fitVonMises([],[],'kappa',simResults{1}.kappa);
for iTuning = 1:length(neuralTuningWidth);
  neuralTuningWidthStr{iTuning} = num2str(neuralTuningWidth(iTuning));
end
figure(f);subplot(2,3,6);
legend(neuralTuningWidthStr);

print('-dpdf',fullfile(figDir,'fig3.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispCategoryFigure    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispCategoryFigure(simResults,dataDir,figDir)

f = mlrSmartfig('testCinvor_category','reuse');clf(f);hold on

subplot(1,4,1);hold on
nCategories = 2;
plot(-90:90,(-90:90>(-180/(2*nCategories))) & (-90:90<(180/(2*nCategories))),'c-','Color',getSmoothColor(1,3,'cool'));
nCategories = 4;
plot(-90:90,(-90:90>(-180/(2*nCategories))) & (-90:90<(180/(2*nCategories))),'c-','Color',getSmoothColor(2,3,'cool'));
nCategories = 8;
plot(-90:90,(-90:90>(-180/(2*nCategories))) & (-90:90<(180/(2*nCategories))),'c-','Color',getSmoothColor(3,3,'cool'));

xaxis(-90,90);
yaxis(0,1);
xlabel('orientation (deg)');
ylabel('Neural response');
drawPublishAxis('xTick',[-90 0 90]);

r2 = 0.2;r2binwidth = 0.1;
subplot(1,4,2);
dispChannelTuning(simResults{1},r2,r2binwidth);
title('exponent = 3');

subplot(1,4,3);
dispChannelTuning(simResults{2},r2,r2binwidth);
title('exponent = 8');

subplot(1,4,4);
dispChannelTuning(simResults{3},r2,r2binwidth)
title('Stick function');

print('-dpdf',fullfile(figDir,'figCategory.pdf'));

% pairedttest
%
%      usage: pairedttest(dist1, dist2)
%       e.g.: p = pairedttest(dist1, dist2)
%         by: ustin gardner
%       date: 2013/02/06
%
function [p t dof d] = pairedttest(dist1,dist2)

% check arguments
if (nargin ~= 2)
  help myttest;
  return;
end

% check length
dist1 = dist1(:);dist2 = dist2(:);
if (length(dist1) ~= length(dist2))
  disp(sprintf('(pairedttest) Distributions must be the same size'));
  p = nan;
  return
end

% remove nan's from the data
goodvals = find(~isnan(dist1(:)) & ~isnan(dist2(:)));
dist1 = dist1(goodvals);
dist2 = dist2(goodvals);

if (length(dist1) <= 2)
  disp(sprintf('(pairedttest) Distributions must both have 2 or more non-nan matching values'));
  p = nan;
  return
end

% number of values
n = length(dist1);

% calculate mean paired difference
md = mean(dist1-dist2);

% calculate variances
var1=var(dist1);
var2=var(dist2);

% standard error of the mean
Sd = std(dist1-dist2)/sqrt(n);

% t statistic
t = md/Sd;

% p-value
dof = n-1;
p = betainc(dof/(dof+t^2),dof/2,.5);

% return some stuff
if nargout == 4
  d.p = p;
  d.n = n;
  d.t = t;
  d.m1 = mean(dist1);
  d.m2 = mean(dist2);
  d.std1 = sqrt(var1);
  d.std2 = sqrt(var2);
end


