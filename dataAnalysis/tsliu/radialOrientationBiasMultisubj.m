subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
subjNameList = {'s005','s007','s017','s019','s021','s029'};
nsubj=length(subjList);
analName='decon1gIns';
left_anal_name = {'pRF_left','pRF-new','pRF_V1','pRF-left','pRF-left','pRF-left'};
right_anal_name = {'pRF_right','pRF-new','pRF_V1','pRF-right','pRF','pRF-right'};
left_anal_var_name = {'pRF_left','pRF_new','pRF_V1','pRF_left','pRF_left','pRF_left'};
ROI = 2;
nsubj=6;
M = 8;
carlsonAnal = zeros(2,2,M); %angularZone/Contrast/Orientation
biasAnal = zeros(2,2,M); %angularZone/Contrast/Orientation
for subj = 1:6
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

	for iValue = 1:2
		importantVoxels = prf_voxels(1,:)>0.1;
		x_cutoff = -8;
		y_cutoff = 0;
		if(iValue == 1)
			importantVoxels = importantVoxels & (prf_voxels(3,:).*cos(prf_voxels(2,:)) - x_cutoff).*(prf_voxels(3,:).*sin(prf_voxels(2,:)) - y_cutoff) < 0;
		else
			importantVoxels = importantVoxels & (prf_voxels(3,:).*cos(prf_voxels(2,:)) - x_cutoff).*(prf_voxels(3,:).*sin(prf_voxels(2,:)) - y_cutoff) > 0;
		end
		orientation = 0:22.5:179;
		N = s.lvf.roi{ROI}.n;
		M = 8;
		voxelResponse = zeros(N,M);
		for i = 1:M
			voxelResponse(:,i) = mean(s.lvf.roi{ROI}.instance.instances{i},1)';
		end
		avgTransform = mean(voxelResponse(importantVoxels,:),1);
		carlsonAnal(iValue,1,:) = carlsonAnal(iValue,1,:) + reshape(1.0/nsubj*avgTransform,[1,1,M]);

		voxelResponse = zeros(N,M);
		for i = 1:M
			voxelResponse(:,i) = mean(s.lvf.roi{ROI}.instance.instances{i+M},1)';
		end
		avgTransform = mean(voxelResponse(importantVoxels,:),1);
		%avgTransform = min(1,avgTransform); %makes plots nicer; 
		carlsonAnal(iValue,2,:) = carlsonAnal(iValue,2,:) + reshape(1.0/nsubj*avgTransform,[1,1,M]);
	end

	for iValue = 1:2
		importantVoxels = prf_voxels(1,:)>0.1;
		x_cutoff = -8;
		y_cutoff = 0;
		if(iValue == 1)
			importantVoxels = importantVoxels & (prf_voxels(3,:).*sin(prf_voxels(2,:)) - y_cutoff) < 0;
		else
			importantVoxels = importantVoxels & (prf_voxels(3,:).*sin(prf_voxels(2,:)) - y_cutoff) > 0;
		end
		orientation = 0:22.5:179;
		N = s.lvf.roi{ROI}.n;
		M = 8;
		voxelResponse = zeros(N,M);
		for i = 1:M
			voxelResponse(:,i) = mean(s.lvf.roi{ROI}.instance.instances{i},1)';
		end
		avgTransform = mean(voxelResponse(importantVoxels,:),1);
		biasAnal(iValue,1,:) = biasAnal(iValue,1,:) + reshape(1.0/nsubj*avgTransform,[1,1,M]);

		voxelResponse = zeros(N,M);
		for i = 1:M
			voxelResponse(:,i) = mean(s.lvf.roi{ROI}.instance.instances{i+M},1)';
		end
		avgTransform = mean(voxelResponse(importantVoxels,:),1);
		%avgTransform = min(1,avgTransform); %makes plots nicer; 
		biasAnal(iValue,2,:) = biasAnal(iValue,2,:) + reshape(1.0/nsubj*avgTransform,[1,1,M]);
	end
end

figure;
for iValue = 1:2
	subplot(1,2,iValue);
	plot(orientation,reshape(carlsonAnal(iValue,1,:),[1,M]),'k');
	hold;
	plot(orientation,reshape(carlsonAnal(iValue,2,:),[1,M]),'r');
end

figure;
for iValue = 1:2
	subplot(1,2,iValue);
	plot(orientation,reshape(biasAnal(iValue,1,:),[1,M]),'k');
	hold;
	plot(orientation,reshape(biasAnal(iValue,2,:),[1,M]),'r');
end
