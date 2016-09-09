% testCinvor.m
%
%      usage: testCinvor()
%         by: justin gardner
%       date: 09/07/16
%    purpose: function to test build/test channels based on cinvor data
%
function retval = testCinvor(varargin)

% get arguments
getArgs(varargin,{'dataDir=~/data/cinvor2','analName=decon1gIns','contrastName=high','nFold=5','subjectList',{'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'}});

% number of subjects
nSubjects = length(subjectList);

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
      leftFoldChannel(iFold) = buildChannels(trainInstances,e.stimVals,'fitNoiseModel=1','noiseModelGridSearchOnly=1');
      leftFoldOutputs(iFold) = testChannels(testInstances,e.stimVals,leftFoldChannel(iFold));
      

      % comupte train/test of right visual field channels
      [trainInstances testInstances] = getCrossValInstances(rvfInstances,crossVal,iFold);
      rightFoldChannel(iFold) = buildChannels(trainInstances,e.stimVals,'fitNoiseModel=1','noiseModelGridSearchOnly=1');
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

% display figure
mlrSmartfig(sprintf('testCinvor%s',contrastName),'reuse');clf;

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

keyboard

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
  mylegend({'Contra visual field','Ipsi visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
end

subplot(2,nCol,iCol+nCol);
statsStrLeft = dispChannelOutput(leftOutput,leftChannel,'likelihood=1');
statsStrRight = dispChannelOutput(rightOutput,rightChannel,'MarkerFaceColor=r','likelihood=1');
title(sprintf('Likelihood\n (left: %s)\n(right: %s)',statsStrLeft,statsStrRight));
xlabel('Orientation difference from true (deg)');ylabel('Likelihood (p)');
if ~contraIpsi
  mylegend({'Right visual field','Left visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
else
  mylegend({'Contra visual field','Ipsi visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
end

drawnow