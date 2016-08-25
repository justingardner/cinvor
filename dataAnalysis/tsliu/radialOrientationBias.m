mrQuit;
subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
subjNameList = {'s005','s007','s017','s019','s021','s029'};
nsubj=length(subjList);
analName='decon1gIns';
left_anal_name = {'pRF_left','pRF-new','pRF_V1','pRF-left','pRF-left','pRF-left'};
right_anal_name = {'pRF_right','pRF-new','pRF_V1','pRF-right','pRF','pRF-right'};
left_anal_var_name = {'pRF_left','pRF_new','pRF_V1','pRF_left','pRF_left','pRF_left'};
right_anal_var_name = {'pRF_right','pRF_new','pRF_V1','pRF_right','pRF','pRF_right'};
ROI = 2;
nsubj=6;
M = 8;
subj = 1;
load([subjList{subj},'Anal/',analName]);
subjName = subjNameList{subj};
%coordinates of voxels
s.lvf.roi{ROI}.scanCoords;
if(ROI == 1)
	myData = load(['/Users/gru/data/ganatret/',subjName,'/retinotopy_',subjName,'/Averages/pRFAnal/',left_anal_name{subj},'.mat']);
	pRF = myData.(left_anal_var_name{subj});
else
	myData = load(['/Users/gru/data/ganatret/',subjName,'/retinotopy_',subjName,'/Averages/pRFAnal/',right_anal_name{subj},'.mat']);
	pRF = myData.(right_anal_var_name{subj});
end
r2_mat = pRF.overlays(1).data{end};
nVoxels = s.lvf.roi{ROI}.n;
transMat = generateTransitionMatrix(['/Users/gru/data/ganatret/',subjName,'/retinotopy_',subjName],['/Users/gru/data/cinvor2/',subjList{subj}],1);
myCoords = [s.lvf.roi{ROI}.scanCoords; ones(1,nVoxels)];
myCoords = transMat*myCoords;
myCoords = round(myCoords(1:3,:));
N_prf = 4;
prf_voxels = zeros(N_prf,nVoxels);
for j = 1:N_prf
	this_mat = pRF.overlays(j).data{end};
	for i = 1:nVoxels
		coords = myCoords(:,i);
		prf_voxels(j,i) = this_mat(coords(1),coords(2),coords(3));
	end
end
figure;
importantVoxels = prf_voxels(1,:)>0.5;
plot(prf_voxels(3,importantVoxels).*cos(prf_voxels(2,importantVoxels)),prf_voxels(3,importantVoxels).*sin(prf_voxels(2,importantVoxels)),'o');
axis([-20 20 -20 20]);

myR2 = zeros(1,nVoxels);
instanceMat = zeros(0,nVoxels);
for i = 1:M
	instanceMat = [instanceMat; s.lvf.roi{ROI}.instance.instances{i+M}];
end
totalVar = var(instanceMat,1);
residualVar = zeros(1,nVoxels);
for i = 1:M
	residualVar = residualVar + 1.0/M*var(s.lvf.roi{ROI}.instance.instances{i+M},1)
end
myR2 = 1 - residualVar./totalVar;

myR2 = zeros(1,nVoxels);
instanceMat = zeros(0,nVoxels);
for i = 1:2*M
	instanceMat = [instanceMat; s.lvf.roi{ROI}.instance.instances{i}];
end
totalVar = var(instanceMat,1);
residualVar = zeros(1,nVoxels);
for i = 1:2*M
	residualVar = residualVar + 1.0/(2*M)*var(s.lvf.roi{ROI}.instance.instances{i},1)
end
myR2 = 1 - residualVar./totalVar;


figure;
for iValue = 1:4
	subplot(2,2,iValue);
	importantVoxels = prf_voxels(1,:)>0.2 & s.lvf.roi{ROI}.r2 > 0.15;
	%plot(prf_voxels(3,importantVoxels).*cos(prf_voxels(2,importantVoxels)),prf_voxels(3,importantVoxels).*sin(prf_voxels(2,importantVoxels)),'o');


	%test carlson idea
	x_cutoff = -8;
	y_cutoff = 0;
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
figure;
for iValue = 1:2
	subplot(1,2,iValue);
	importantVoxels = prf_voxels(1,:)>0.1;
	prf_voxels(2,importantVoxels) = 2*pi*(prf_voxels(2,importantVoxels) < 0) + prf_voxels(2,importantVoxels);
	cutoff = median(prf_voxels(2,importantVoxels)); %median
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