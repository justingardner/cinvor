% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=16','trialPerStim=24');


scenarioName = 'Tuning width modulation only';
testParameter = 'kappa';
testParameterValues = [10];


% display which scenario we are running
disp(scenarioName);
N_trials = 24;
N = 100;
M = 16;
m = setCinvorModel(e,testParameter,testParameterValues(1));
trainInstances = getCinvorInstances(m,e);
if(false)
	trainMatrix = zeros(N_trials*M,N);
	for i = 1:M
		trainMatrix(1 + (i-1)*N_trials:i*N_trials,:) = trainInstances{i};
	end
	[trainMatrix, mu, sigma] = zscore(trainMatrix);
	for i = 1:M
		trainInstances{i} = trainMatrix(1 + (i-1)*N_trials:i*N_trials,:);
	end
end

voxelResponse = zeros(N,M);
for i = 1:M
	voxelResponse(:,i) = mean(trainInstances{i},1)';
end
%voxelResponse = zscore(voxelResponse,0,2); %zScore
F = fft(voxelResponse,[],2); %takes Fourier transform of each row.
%
%normConst = mean(F(:,1));
%F = F./normConst; %normalize by first component

%fit noise model
std_mat = zeros(M,N);
for i = 1:M
	std_mat(i,:) = std(trainInstances{i});%/normConst
end
thisVar = (mean(std_mat)/sqrt(N_trials)).^2;
FourierVariance = thisVar*M;

observedMean = var(F,[],1);
observedVar1 = N/(N-1)*moment(real(F),2,1);
observedVar2 = N/(N-1)*moment(imag(F),2,1);
observedVariance = observedVar1 + observedVar2;
observedFourth1 = (N^3*moment(real(F),4,1) - (N-1)*3*(2*N-3)*observedVar1.^2)/((N-1)*(N^2-3*N+3)); %fourthmoment
observedFourth2 = (N^3*moment(imag(F),4,1) - (N-1)*3*(2*N-3)*observedVar2.^2)/((N-1)*(N^2-3*N+3)); %fourthmoment
observedVarofVar1 = observedFourth1/N - observedVar1.^2/N*(N-3)/(N-1);
observedVarofVar2 = observedFourth2/N - observedVar2.^2/N*(N-3)/(N-1);
totalVar = fftshift(observedVarofVar1 + observedVarofVar2); %TODO: add VarofFourierVariance
%clf

FourierSize = abs(F).^2 - repmat(FourierVariance',1,M);
FourierSize(:,1) = 0;
%FourierSize = FourierSize./repmat(FourierSize(:,2),1,M);
plotVect = fftshift(mean(FourierSize));
%plot(plotVect)
%axis([1 8 -0.5e-3 3e-3]);


%avgTransform = fftshift(mean(abs(F),1));
%avgTransform = min(1,avgTransform); %makes plots nicer; 
%clf;
%plot(avgTransform); 