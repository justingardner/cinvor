global myRealChannel taoParam sigmaParam rhoParam
myRealChannel = realChannel;
% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=27');

correctionFactor = 4;
M = 8;
trainInstances = cell(1,M);
nVoxels = size(myRealChannel.channelWeights,2);
N_trials = 24;
mean_response = myRealChannel.channelWeights'*myRealChannel.idealStimResponse;
rho = myRealChannel.rho;
tao = 2*myRealChannel.tao;
sigma = 2*myRealChannel.sigma;
taoParam = tao;
sigmaParam = sigma;
rhoParam = rho;
thisError = zeros(M*N_trials,nVoxels);
omega = rho*tao'*tao + (1-rho)*diag(diag(tao'*tao))+sigma^2*myRealChannel.channelWeights'*myRealChannel.channelWeights; 
for i = 1:M
  %trainInstances{i} = mvnrnd(mean_response(:,i),correctionFactor*myRealChannel.omega,N_trials);
  trainInstances{i} = mvnrnd(mean_response(:,i),omega,N_trials);
  thisError(1+(i-1)*N_trials:i*N_trials,:) = trainInstances{i} - repmat(mean_response(:,i),1,N_trials)';
end
thisError2 = mvnrnd(zeros(nVoxels,1),omega,N_trials*M);
thisError3 = mvnrnd(zeros(nVoxels,1),omega,N_trials*M);
diag(1/(N_trials*M)*(thisError3'*thisError3 - thisError2'*thisError2))
diag(1/(N_trials*M)*(thisError'*thisError - thisError2'*thisError2))
diag(1/(N_trials*M)*(thisError'*thisError - trueError'*trueError))

channel = buildChannelsTesting(trainInstances,e.stimVals);


