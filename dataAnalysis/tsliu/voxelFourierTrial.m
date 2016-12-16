
subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
nsubj=length(subjList);
analName='decon1gIns';
load([subjList{5},'Anal/',analName]);
s.lvf.roi{1}.instance.instances;
orientation = 0:22.5:179;
N = s.lvf.roi{1}.n;
M = 8;
N_trials = size(s.lvf.roi{1}.instance.instances{1},1);
voxelResponse = zeros(N,M);
for i = 1:M
	voxelResponse(:,i) = mean(s.lvf.roi{1}.instance.instances{i},1)';
end
F = fft(voxelResponse,[],2); %takes Fourier transform of each row.

%fit noise model
std_mat = zeros(M,N);
for i = 1:M
	std_mat(i,:) = std(s.lvf.roi{1}.instance.instances{i});%/normConst
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
avgSTD = sqrt(totalVar);
%avgTransform = min(1,avgTransform); %makes plots nicer; 
clf
plotVect(5) = 0;
avgSTD(5) = 0;
plot(plotVect);
hold on;
errorbar(1:M,plotVect,avgSTD*1.96)
%axis([0 M+1 -0.5e-3 3e-3]);
