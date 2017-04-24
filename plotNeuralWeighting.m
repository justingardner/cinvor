% simuCinvor.m
%
%      usage: simuTuningChange(scenario)
%         by: dylan cable
%       date: 7/2016
%    purpose: Simulation of the forward encoding model
%             in which model is trained on one tuning 
%             width and tested on another:
%
%             1 = amplitude simulation
%             2 = noise simulation
%             3 = tuning width simulation
%
function FWHM = simuTuningChange(firstKappa,varargin)

%

% get arguments
getArgs(varargin,'dispFig=0');

% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=24');


testParameter = 'kappa';


% display which scenario we are running

m = setCinvorModel(e,testParameter,firstKappa);

trainInstances = getCinvorInstances(m,e);
instanceMatrix=[];
for istim=1:length(trainInstances)
  instanceMatrix=[instanceMatrix; trainInstances{istim}];
end
roundedStimVal = round(e.allStimVal-0.1);
stimMatrix = zeros(180,192);
for i = 1:192
  stimMatrix(roundedStimVal(i)+1,i) = 1;
end
instances = m.ws*m.neuralResponse*stimMatrix;

channel=buildChannels(trainInstances,e.stimVals,'dispChannels=0');
channelNeuralWeights = (channel.channelWeights*m.ws);
channelResponse = channelNeuralWeights*m.neuralResponse;
