% testCinvor.m
%
%      usage: testCinvor()
%         by: justin gardner
%       date: 09/07/16
%    purpose: function to test build/test channels based on cinvor data
%
function retval = testCinvor(varargin)

% check arguments
if ~any(nargin == [0])
  help testCinvor
  return
end

% get arguments
getArgs(varargin,{'dataDir=~/data/cinvor','analName=decon1gIns','contrastName=high','nFold=5','subjectList',{'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'}});

% number of subjects
nSubjects = length(subjectList);

% display what we are doing
disppercent(-inf,sprintf('(testCinvor) Building and testing channels for %s contrast using analysis %s with %i fold cross-validation on %i subjects',contrastName,analName,nFold,nSubjects));

% iterate over subjects
for iSubject = 1:nSubjects
  % load the data
  dataName = [subjectList{iSubject},'Anal/',analName];
  disp(sprintf('(testCinvor) Loading data for subject %i/%i: %s',iSubject,nSubjects,dataName));
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
      channel = buildChannels(trainInstances,e.stimVals);
      leftFoldOutputs(iFold) = testChannels(testInstances,e.stimVals,channel);
      %leftOutputs(iSubject,iROI,iFold) = testChannels(testInstances,e.stimVals,channel);
      % comupte train/test of right visual field channels
      [trainInstances testInstances] = getCrossValInstances(rvfInstances,crossVal,iFold);
      channel = buildChannels(trainInstances,e.stimVals);
      rightFoldOutputs(iFold) = testChannels(testInstances,e.stimVals,channel);
      %rightOutputs(iSubject,iROI,iFold) = testChannels(testInstances,e.stimVals,channel);
      disppercent(calcPercentDone(iSubject,nSubjects,iROI,nROI,iFold,nFold));
    end
    
    % average over folds and store in structure
    leftOutputs(iSubject,iROI) = combineChannelOutputs(leftFoldOutputs,'combineAveraged');
    rightOutputs(iSubject,iROI) = combineChannelOutputs(rightFoldOutputs,'combineAveraged');
  end
end
disppercent(inf);


% display figure
mlrSmartfig(sprintf('testCinvor%s',contrastName),'reuse');clf;

for iROI = 1:2:nROI
  % get contras/ipsi for left
  contraLeft = combineChannelOutputs(rightOutputs(:,iROI));
  ipsiLeft = combineChannelOutputs(leftOutputs(:,iROI));

  % get contras/ipsi for right
  contraRight = combineChannelOutputs(leftOutputs(:,iROI+1));
  ipsiRight = combineChannelOutputs(rightOutputs(:,iROI+1));

  % get combined contra / ipsi
  contraCombined = combineChannelOutputs([rightOutputs(:,iROI) leftOutputs(:,iROI+1)],channel);
  ipsiCombined = combineChannelOutputs([rightOutputs(:,iROI+1) leftOutputs(:,iROI)],channel);

  % plot left ROI
  subplot(1,nROI+nROI/2,iROI);
  dispChannelOutput(contraLeft,channel);
  dispChannelOutput(ipsiLeft,channel,'MarkerFaceColor=r');

  % set title and legends
  title(sprintf('%s contrast: %s (%i fold cross-validation, nSubjects = %i)',roiName{iROI},contrastName,nFold,nSubjects));
  mylegend({'Left visual field','Right visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
  ylabel('Channel response (percentile of full response');

  % plot right ROI
  subplot(1,nROI+nROI/2,iROI+1);
  dispChannelOutput(contraRight,channel,'MarkerFaceColor=r');
  dispChannelOutput(ipsiRight,channel);

  % set title and legends
  title(sprintf('%s contrast: %s (%i fold cross-validation, nSubjects = %i)',roiName{iROI+1},contrastName,nFold,nSubjects));
  mylegend({'Left visual field','Right visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
  xlabel('Orientation difference from actual (degrees)');

  % plot combined ROI
  subplot(1,nROI+nROI/2,iROI+2);
  dispChannelOutput(contraCombined,channel,'MarkerFaceColor=r');
  dispChannelOutput(ipsiCombined,channel);

  % set title and legends
  title(sprintf('%s contrast: %s (%i fold cross-validation, nSubjects = %i)',roiName{iROI+1}(2:end),contrastName,nFold,nSubjects));
  mylegend({'Contra visual field','Ipsi visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});

  makeEqualYaxis(1,nROI+nROI/2,iROI:iROI+2);
end

keyboard

