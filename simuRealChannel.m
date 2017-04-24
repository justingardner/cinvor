global myRealChannel taoParam sigmaParam rhoParam
myRealChannel = realChannel;

% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=27');
T = 1;
totalTestResponse = zeros(1,9);
correctionFactor = 4;
rho = 0.9*myRealChannel.rho;
tao = 0.95*myRealChannel.tao;
sigma = 1.5*myRealChannel.sigma;
omega = rho*tao'*tao + (1-rho)*diag(diag(tao'*tao))+sigma^2*myRealChannel.channelWeights'*myRealChannel.channelWeights; 

taoParam = tao;
sigmaParam = sigma;
rhoParam = rho;
M = 8;
trainInstances = cell(1,M);
nVoxels = size(myRealChannel.channelWeights,2);
N_trials = 24;%21
mean_response = myRealChannel.channelWeights'*myRealChannel.idealStimResponse;
for i = 1:M
  %trainInstances{i} = mvnrnd(mean(realTrainInstances{i})',myRealChannel.omega,N_trials);
  trainInstances{i} = mvnrnd(mean_response(:,i),omega,N_trials);
end
N_trials_test = 6;
testInstances = cell(1,M);
contrast2 = 0.25;
for i = 1:M
  %testInstances{i} = mvnrnd(mean(realTestInstances{i})',myRealChannel.omega,N_trials_test);
  testInstances{i} = mvnrnd(contrast2*mean_response(:,i),omega,N_trials_test);
end

channel = buildChannels(trainInstances,e.stimVals);
%channel = buildChannelsTesting(trainInstances,e.stimVals);

[avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel);

avgTestResponse = [avgTestResponse, avgTestResponse(1)];
totalTestResponse = totalTestResponse + 1/T*avgTestResponse;

clf;
plot(totalTestResponse);

%strength of signal is about the same **
corr(mean(trainInstances{1})',mean(testInstances{1})')
corr(mean(realTrainInstances{1})',mean(realTestInstances{1})')
%compare std
mean(std(trainInstances{2}))
mean(std(realTrainInstances{2}))
%compate channel noise params
mean(diag(channel.omega))
mean(diag(myRealChannel.omega))

mean(diag(channel.sigma))
mean(diag(myRealChannel.sigma))

%difference in channelWeights! std but not mean
mean(std(channel.channelWeights))
mean(std(myRealChannel.channelWeights))

mean(mean(channel.channelWeights))
mean(mean(myRealChannel.channelWeights))