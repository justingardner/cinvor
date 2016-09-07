% testCinvor.m
%
%      usage: testCinvor()
%         by: justin gardner
%       date: 09/07/16
%    purpose: function to test build/test channels based on cinvor data
%
function retval = testCinvor()

% check arguments
if ~any(nargin == [0])
  help testCinvor
  return
end

% data directories
dataDir = '~/data/cinvor';
subjectList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};

% number of subjects
nSubjects = length(subjectList);

% which analysis type to use
analName = 'decon1gIns';

% contrast to run (either high or low)
contrast = 'high';

% how many fold cross-validation to use
nFold = 10;

for iSubject = 1:nSubjects
  % load the data
  dataName = [subjectList{iSubject},'Anal/',analName];
  disp(sprintf('(testCinvor) Loading data: %s',dataName));
  load(dataName);
  
  % specify the experiment
  e = specifyCinvorExperiment('stimLevel=8','trialPerStim=1');
 
  % get number of ROIS
  nROI = length(s.lvf.roi);

  % decide on which contrast to do
  if strcmp(lower(contrast),'high');
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
    crossVal = getCrossValSets(lvfInstances);
      
    % now train and test over folds
    for iFold = 1:5
      % comupte train/test of left visual field channels
      [trainInstances testInstances] = getCrossValInstances(lvfInstances,crossVal,iFold);
      channel = buildChannels(trainInstances,e.stimVals);
      leftChannelResponse(iSubject,iROI,iFold,:) = testChannels(testInstances,e.stimVals,channel);
      % comupte train/test of right visual field channels
      [trainInstances testInstances] = getCrossValInstances(rvfInstances,crossVal,iFold);
      channel = buildChannels(trainInstances,e.stimVals);
      rightChannelResponse(iSubject,iROI,iFold,:) = testChannels(testInstances,e.stimVals,channel);
    end
  end
end

% display output
mlrSmartfig('testCinvor','reuse');clf;
xVals = circshift(e.stimVals(:),4);
xVals(1:4) = xVals(1:4)-180;
for iROI = 1:nROI
  % select subplot
  subplot(1,nROI,iROI);
  % average across subjects and fold
  thisLeftChannelResponse = leftChannelResponse(:,iROI,:,:);
  thisLeftChannelResponseSTE = squeeze(std(squeeze(mean(thisLeftChannelResponse,3)))/sqrt(nSubjects));
  thisLeftChannelResponse = squeeze(mean(mean(thisLeftChannelResponse,1),3));
  thisRightChannelResponse = rightChannelResponse(:,iROI,:,:);
  thisRightChannelResponseSTE = squeeze(std(squeeze(mean(thisRightChannelResponse,3)))/sqrt(nSubjects));
  thisRightChannelResponse = squeeze(mean(mean(thisRightChannelResponse,1),3));
  % display
  myerrorbar(xVals,thisLeftChannelResponse,'yError',thisLeftChannelResponseSTE,'MarkerFaceColor','r','MarkerEdgeColor','w','MarkerSize',8);hold on
  myerrorbar(xVals,thisRightChannelResponse,'yError',thisRightChannelResponseSTE,'MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',8);hold on
  % set title
  title(sprintf('%s contrast: %s (%i fold cross-validation)',roiName{iROI},contrast,nFold));
  % set legend
  mylegend({'Left visual field','Right visual field'},{{'ro' 'MarkerFaceColor=r' 'MarkerEdgeColor=w'},{'ko' 'MarkerFaceColor=k' 'MarkerEdgeColor=w'}});
  xlabel('Orientation difference from actual (degrees)');
  ylabel('Channel response (a.u.)');
end
  
  


keyboard
