subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
nsubj=length(subjList);
analName='decon1gIns';

nROI = 2;
nSubj = 6;
nCond = 2;
totalMeans = zeros(nROI,nSubj,nCond);
%compute means
for subj = 1:nSubj
	load([subjList{subj},'Anal/',analName]);
	for ROI = 1:nROI
		nVoxels = s.rvf.roi{ROI}.n;
		trainInstances = s.lvf.roi{ROI}.instance.instances(1:M);
		for i = 1:M
			totalMeans(ROI,subj,1) = totalMeans(ROI,subj,1) + 1/M*mean(mean(trainInstances{i}));
		end

		trainInstances = s.lvf.roi{ROI}.instance.instances(M+1:2*M);
		for i = 1:M
			totalMeans(ROI,subj,2) = totalMeans(ROI,subj,2) + 1/M*mean(mean(trainInstances{i}));
		end
	end
end

%correlation between low contrast and high contrast
ROI = 1;
nVoxels = 78;
M = 8;
N_trials = 27;
stopN = 27;
nROI = 2;
nSubj = 6;
totalCorr = zeros(nROI,nSubj,N_trials);
Trials = 1;
for ROI = 1:nROI
	for subj = 1:nSubj
		load([subjList{subj},'Anal/',analName]);
		nVoxels = s.rvf.roi{ROI}.n;
		for stopN = N_trials:N_trials
			for T = 1:Trials
				voxelMeans = zeros(M,nVoxels);
				voxelStderr = zeros(M,nVoxels);
				trainInstances = s.rvf.roi{ROI}.instance.instances(1:M);
				for i = 1:M
					voxelMeans(i,:) = mean(trainInstances{i}(randsample(N_trials,stopN),:));
					voxelStderr(i,:) = std(trainInstances{i}(randsample(N_trials,stopN),:))/sqrt(N_trials);
				end


				trainInstances = s.rvf.roi{ROI}.instance.instances(M+1:2*M);
				M = 8;
				N_trials = 27;
				voxelMeans2 = zeros(M,nVoxels);
				voxelStderr2 = zeros(M,nVoxels);
				for i = 1:M
					voxelMeans2(i,:) = mean(trainInstances{i}(randsample(N_trials,stopN),:));
					voxelStderr2(i,:) = std(trainInstances{i}(randsample(N_trials,stopN),:))/sqrt(N_trials);
				end

				
				voxelToPlot = 1;
				clf;
				plot(1:M,voxelMeans2(:,voxelToPlot),'o');
				hold on
				errorbar(1:M,voxelMeans2(:,voxelToPlot),1.96*voxelStderr(:,voxelToPlot),'o');

				totalCorr(ROI,subj,stopN) = totalCorr(ROI,subj,stopN) +  1.0/Trials*mean(diag(corr(voxelMeans,voxelMeans2)));
			end
		end
	end
end
figure;
plot(totalCorr)
xlabel('number of trials');
ylabel('correlation between mean voxel tuning function of low and high contrast');


%correlation between two halves of data in one condition.
ROI = 1;
nVoxels = 78;
M = 8;
stopN = 27;
nROI = 2;
nSubj = 6;
nCond = 2; %high / low 
totalCorr = zeros(nROI,nSubj,nCond);
Trials = 1;
for ROI = 1:nROI
	for subj = 1:nSubj
		load([subjList{subj},'Anal/',analName]);
		nVoxels = s.lvf.roi{ROI}.n;
		for condi = 1:nCond
			voxelMeans = zeros(M,nVoxels);
			voxelMeans2 = zeros(M,nVoxels);
			if(condi == 1)
				trainInstances = s.lvf.roi{ROI}.instance.instances(1:M);
			else
				trainInstances = s.lvf.roi{ROI}.instance.instances(M+1:2*M);
			end
			for i = 1:M
				voxelMeans(i,:) = mean(trainInstances{i}(1:round(N_trials/2),:));
				voxelMeans2(i,:) = mean(trainInstances{i}(round(N_trials/2)+1:end,:));
			end
			totalCorr(ROI,subj,condi) = mean(diag(corr(voxelMeans,voxelMeans2)));
		end
	end
end

varEstim = zeros(nROI,nSubj);
for ROI = 1:nROI
	for subj = 1:nSubj
		load([subjList{subj},'Anal/',analName]);
		trainInstances = s.lvf.roi{ROI}.instance.instances(1+M:M+M);
		nVoxels = s.lvf.roi{ROI}.n;;
		M = 8;
		N_trials = 27;
		voxelMeans = zeros(M,nVoxels);
		voxelStderr = zeros(M,nVoxels);
		for i = 1:M
			voxelMeans(i,:) = mean(trainInstances{i});
			voxelStderr(i,:) = std(trainInstances{i})/sqrt(stopN);
		end
		%subtract out mean
		voxelMeansZ = voxelMeans - repmat(mean(voxelMeans),M,1);%mean(mean(voxelMeans));%repmat(mean(voxelMeans),M,1);
		varEstim(ROI,subj) = mean(mean(8/7*voxelMeansZ.^2 - voxelStderr.^2,2));
	end
end
base_means = mean(voxelMeans);



T = 1000;
all_counts = zeros(T,1);
%all_mean_zero
for count = 1:T
	trueStd = 0.10;
	noiseStd = 0.217*sqrt(N_trials);
	trueMeans = randn(M,nVoxels)*trueStd ;
	genInstances = cell(1,M);
	for i = 1:M
		genInstances{i} = randn(size(trainInstances{i}))*noiseStd + repmat(trueMeans(i,:),N_trials,1);
	end

	voxelMeans = zeros(M,nVoxels);
	voxelStderr = zeros(M,nVoxels);
	for i = 1:M
		voxelMeans(i,:) = mean(genInstances{i});
		voxelStderr(i,:) = std(genInstances{i})/sqrt(N_trials);
	end
	%subtract out mean
	%voxelMeansZ = voxelMeans - repmat(mean(voxelMeans),M,1);%mean(mean(voxelMeans));%repmat(mean(voxelMeans),M,1);%mean(mean(voxelMeans));%repmat(mean(voxelMeans),M,1);

	%Unbiased estimator for the true var:
	%notice the correction
	all_counts(count) = mean(mean(voxelMeans.^2 - voxelStderr.^2,2));
end
sqrt(mean(all_counts))


T = 1000;
all_counts = zeros(T,1);
%different base mean
for count = 1:T
	trueStd = 0.10;
	noiseStd = 0.217*sqrt(N_trials);
	trueMeans = randn(M,nVoxels)*trueStd + repmat(base_means,8,1);
	genInstances = cell(1,M);
	for i = 1:M
		genInstances{i} = randn(size(trainInstances{i}))*noiseStd + repmat(trueMeans(i,:),N_trials,1);
	end

	voxelMeans = zeros(M,nVoxels);
	voxelStderr = zeros(M,nVoxels);
	for i = 1:M
		voxelMeans(i,:) = mean(genInstances{i});
		voxelStderr(i,:) = std(genInstances{i})/sqrt(N_trials);
	end
	voxelMeans = voxelMeans - repmat(mean(voxelMeans),M,1);
	%subtract out mean
	%voxelMeansZ = voxelMeans - repmat(mean(voxelMeans),M,1);%mean(mean(voxelMeans));%repmat(mean(voxelMeans),M,1);%mean(mean(voxelMeans));%repmat(mean(voxelMeans),M,1);

	%Unbiased estimator for the true var:
	%notice the correction 8/7
	all_counts(count) = mean(mean(8/7*voxelMeans.^2 - voxelStderr.^2,2));
end
sqrt(mean(all_counts))








