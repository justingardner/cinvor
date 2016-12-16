subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
nsubj=length(subjList);
analName='decon1gIns';
ROI = 2;
nsubj=6;
M = 8;
subj = 2;
load([subjList{subj},'Anal/',analName]);
subjName = 's007';
%coordinates of voxels
s.lvf.roi{ROI}.scanCoords;
load(['/Users/gru/data/ganatret/',subjName,'/Averages/pRFAnal/pRF-new.mat']);
r2_mat = pRF_right.overlays(1).data{3};
nVoxels = s.lvf.roi{ROI}.n;

N_prf = 3;
prf_voxels = zeros(N_prf,nVoxels);
for j = 1:N_prf
	this_mat = pRF_right.overlays(j).data{3};
	for i = 1:nVoxels
		coords = s.lvf.roi{ROI}.scanCoords(:,i);
		prf_voxels(j,i) = this_mat(coords(1),coords(2),coords(3));
	end
end
figure;
importantVoxels = prf_voxels(1,:)>0.1;
plot(prf_voxels(3,importantVoxels).*cos(prf_voxels(2,importantVoxels)),prf_voxels(3,importantVoxels).*sin(prf_voxels(2,importantVoxels)),'o');

for iValue = 1:4
	subplot(2,2,iValue);
	importantVoxels = prf_voxels(1,:)>0.1;
	%plot(prf_voxels(3,importantVoxels).*cos(prf_voxels(2,importantVoxels)),prf_voxels(3,importantVoxels).*sin(prf_voxels(2,importantVoxels)),'o');


	%test carlson idea
	x_cutoff = -9.5;
	y_cutoff = 6;
	if(iValue == 1)
		importantVoxels = importantVoxels & (prf_voxels(3,:).*cos(prf_voxels(2,:)) < x_cutoff) & (prf_voxels(3,:).*sin(prf_voxels(2,:)) > y_cutoff);
	elseif(iValue == 2)
		importantVoxels = importantVoxels & (prf_voxels(3,:).*cos(prf_voxels(2,:)) > x_cutoff) & (prf_voxels(3,:).*sin(prf_voxels(2,:)) > y_cutoff);
	elseif(iValue == 3)
		importantVoxels = importantVoxels & (prf_voxels(3,:).*cos(prf_voxels(2,:)) < x_cutoff) & (prf_voxels(3,:).*sin(prf_voxels(2,:)) < y_cutoff);
	else
		importantVoxels = importantVoxels & (prf_voxels(3,:).*cos(prf_voxels(2,:)) > x_cutoff) & (prf_voxels(3,:).*sin(prf_voxels(2,:)) < y_cutoff);
	end
	orientation = 0:22.5:179;
	N = s.lvf.roi{ROI}.n;
	M = 8;
	voxelResponse = zeros(N,M);
	for i = 1:M
		voxelResponse(:,i) = mean(s.lvf.roi{ROI}.instance.instances{i},1)';
	end
	avgTransform = mean(voxelResponse(importantVoxels,:),1);
	plot(orientation,avgTransform,'b');
	hold on;
	voxelResponse = zeros(N,M);
	for i = 1:M
		voxelResponse(:,i) = mean(s.lvf.roi{ROI}.instance.instances{i+8},1)';
	end
	avgTransform = mean(voxelResponse(importantVoxels,:),1);
	%avgTransform = min(1,avgTransform); %makes plots nicer; 

	plot(orientation,avgTransform,'r');
end



%test radialOrientationBias
prf_voxels(1,18) = 0.05; %exclude this voxel.
figure;
for iValue = 1:2
	subplot(1,2,iValue);
	importantVoxels = prf_voxels(1,:)>0.1;
	prf_voxels(2,importantVoxels) = 2*pi*(prf_voxels(2,importantVoxels) < 0) + prf_voxels(2,importantVoxels);
	cutoff = 2.5697; %median
	if(iValue == 1)
		importantVoxels = importantVoxels & (prf_voxels(2,:) < cutoff);
	else
		importantVoxels = importantVoxels & (prf_voxels(2,:) > cutoff);
	end

	orientation = 0:22.5:179;
	N = s.lvf.roi{ROI}.n;
	M = 8;
	voxelResponse = zeros(N,M);
	for i = 1:M
		voxelResponse(:,i) = mean(s.lvf.roi{ROI}.instance.instances{i},1)';
	end
	avgTransform = mean(voxelResponse(importantVoxels,:),1);
	plot(orientation,avgTransform,'b');
	hold on;
	voxelResponse = zeros(N,M);
	for i = 1:M
		voxelResponse(:,i) = mean(s.lvf.roi{ROI}.instance.instances{i+8},1)';
	end
	avgTransform = mean(voxelResponse(importantVoxels,:),1);
	%avgTransform = min(1,avgTransform); %makes plots nicer; 

	plot(orientation,avgTransform,'r');
end