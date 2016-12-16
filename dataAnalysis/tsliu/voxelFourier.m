
subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
nsubj=length(subjList);
analName='decon1gIns';
load([subjList{4},'Anal/',analName]);
s.lvf.roi{1}.instance.instances;
orientation = 0:22.5:179;
N = s.lvf.roi{1}.n;
M = 8;
voxelResponse = zeros(N,M);
for i = 1:M
	voxelResponse(:,i) = mean(s.lvf.roi{1}.instance.instances{i},1)';
end
F = fft(voxelResponse,[],2); %takes Fourier transform of each row.
avgTransform = fftshift(mean(abs(F),1));
%avgTransform = min(1,avgTransform); %makes plots nicer; 
clf;
plot(avgTransform,'b');
hold on;
voxelResponse = zeros(N,M);
for i = 1:M
	voxelResponse(:,i) = mean(s.lvf.roi{1}.instance.instances{i+8},1)';
end
F = fft(voxelResponse,[],2); %takes Fourier transform of each row.
avgTransform = fftshift(mean(abs(F),1));
%avgTransform = min(1,avgTransform); %makes plots nicer; 

plot(avgTransform,'r');

%%compute noise for each condition