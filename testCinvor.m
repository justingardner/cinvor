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
highFit = fitVonMises(xVals,highContrastResponse);
lowFit = fitVonMises(xVals,lowContrastResponse);

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

