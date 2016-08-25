% specifyCinvorExperiment.m
%
%      usage: specifyCinvorExperiment()
%         by: Taosheng Liu
%    purpose: 
%
% Set up an experiment, with a certain number of stimulus values and number
% of trials per stimulus value, default to 8 and 24, respectively.
function e=specifyCinvorExperiment(varargin)

% parse arguments
getArgs(varargin, {'stimLevel=8','trialPerStim=24','plotModel=0','figTitle','Default cinvor model'});

% full range of degrees
e.totalRangeDeg = 180;

% set the number of stimulus levels
e.stimLevel=stimLevel;
e.stimVals=0:e.totalRangeDeg/stimLevel:(e.totalRangeDeg-1);

% number of trials per stimulus
e.trialPerStim=trialPerStim;

% create the array of all stimulus values
allStimVal=[];
for i=1:stimLevel
  allStimVal=[allStimVal, repmat(e.stimVals(i),1,trialPerStim)];
end
%contains stimulus values on each trial
e.allStimVal=allStimVal;

% set whether to plot models
e.plotModel = plotModel;

e.figTitle = figTitle;