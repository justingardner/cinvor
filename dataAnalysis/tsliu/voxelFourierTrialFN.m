
subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
nsubj=length(subjList);
analName='decon1gIns';
plotVectList = zeros(nsubj,8);
for subj = 1:nsubj
	load([subjList{subj},'Anal/',analName]);
	s.lvf.roi{1}.instance.instances;
	orientation = 0:22.5:179;
	N = s.lvf.roi{2}.n;
	M = 8;
	N_trials = size(s.lvf.roi{2}.instance.instances{1},1);
	voxelResponse = zeros(N,M,N_trials);
	for i = 1:M
		voxelResponse(:,i,:) = s.lvf.roi{2}.instance.instances{i}';
	end
	%voxelResponse = zscore(voxelResponse,0,2); %zScore
	F = fft(voxelResponse,[],2); %takes Fourier transform of each row.
	F = real(F);
	%F = cat(1,rand(thisSize),rand(thisSize)+2);
	%F = rand(thisSize);%%+i*rand(thisSize);
	%fit noise model
	std_mat = std(F,0,3);


	%clf
	FourierSize = mean(abs(F).^2,3) - std_mat.^2; %switch mean and abs
	%FourierSize(:,1) = 0;
	%FourierSize = FourierSize./repmat(FourierSize(:,2),1,M);
	plotVect = fftshift(mean(FourierSize));
	plotVectList(subj,:) = plotVect;
end
avgSTD = sqrt(totalVar);
%avgTransform = min(1,avgTransform); %makes plots nicer; 
clf
plotVect(5) = 0;
%avgSTD(5) = 0;
plot(plotVect);
hold on;
%errorbar(1:M,plotVect,avgSTD*1.96)
%axis([0 M+1 -0.5e-3 3e-3]);

G = F(:,2:end,:);
plot(reshape((abs(mean(G,3)).^2),[size(F,1)*7,1]) , reshape(var(G,0,3)/27,[size(F,1)*7,1]),'o');
hold on;
plot(0:8,0:8);


G = 8*randn(128,7,27);
plot(reshape((abs(mean(G,3)).^2),[128*7,1]) , reshape(var(G,0,3)/27,[128*7,1]),'o');
hold on;
plot(0:8,0:8);



F = fft(voxelResponse,[],2); %takes Fourier transform of each row.
F = real(F);
G = F(:,2:end,:);
H = reshape(G,[size(F,1)*7,27])
mean(H,1) %very wierd property in data;: mean across different trials is very different.


%another thing
F = fft(voxelResponse,[],2); %takes Fourier transform of each row.
F = real(F);
G = F(:,2:end,:);
std_mat = std(G,0,3);
FourierSize = mean(abs(G).^2,3) - std_mat.^2; %switch mean and abs
mean(mean(FourierSize)) %negative
H = reshape(G,[size(F,1)*7,27]);

mean(mean(H.^2)) - mean(std(H,0,2))^2; %positive
mean(mean(H.^2)) - mean(std(H,0,2).^2); %negative.  In fact, equal to mean(mean(FourierSize))
