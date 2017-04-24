

% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=27');
T = 100;
totalTestResponse = zeros(1,9);
correctionFactor = 0.3;
for trial = 1:T
	rho = 0*myRealChannel.rho;
	tao = 21/20*myRealChannel.tao;
	sigma = 0*1.5*myRealChannel.sigma;
	omega = rho*tao'*tao + (1-rho)*diag(diag(tao'*tao))+sigma^2*myRealChannel.channelWeights'*myRealChannel.channelWeights; 

	taoParam = tao;
	sigmaParam = sigma;
	rhoParam = rho;
	M = 8;
	trainInstances = cell(1,M);
	nVoxels = size(myRealChannel.channelWeights,2);
	N_trials = 21;%21
	mean_response = correctionFactor*myRealChannel.channelWeights'*myRealChannel.idealStimResponse;
	for i = 1:M
	  %trainInstances{i} = mvnrnd(mean(realTrainInstances{i})',myRealChannel.omega,N_trials);
	  trainInstances{i} = mvnrnd(mean_response(:,i),omega,N_trials);
	end
	N_trials_test = 5;
	testInstances = cell(1,M);
	contrast2 = 1;%1/1.74;
	for i = 1:M
	  %testInstances{i} = mvnrnd(mean(realTestInstances{i})',myRealChannel.omega,N_trials_test);
	  testInstances{i} = mvnrnd(contrast2*correctionFactor*mean_response(:,i),omega,N_trials_test);
	end

	channel = buildChannels(trainInstances,e.stimVals,'fitNoise',0);
	%channel = buildChannelsTesting(trainInstances,e.stimVals);

	[avgTestResponse r2 classifyCorrTotal stimValVector predStimVal posterior] = testChannels(testInstances,e.stimVals,channel,'fitNoise',0);
	avgTestResponse = [avgTestResponse, avgTestResponse(1)];
	totalTestResponse = totalTestResponse + 1/T*avgTestResponse;
end
clf;
plot(totalTestResponse);


mean(std(channel.channelWeights))
mean(std(myRealChannel.channelWeights))

corr(channel.channelWeights(1,:)',channel.channelWeights(5,:)')
corr(myRealChannel.channelWeights(1,:)',myRealChannel.channelWeights(5,:)')