all_corrs = zeros(1,S);
S = 1000;
for i = 1:S
	instances = randn(26,8);
	mean_corr = 0;
	T = 1;
	for j = 1:T
		shuffleInstances = instances(randperm(26),:);
		mean_corr = mean_corr + 1/T*corr(mean(shuffleInstances(1:13,:))',mean(shuffleInstances(14:26,:))');
	end
	all_corrs(i) = mean_corr;
end
sorted_corrs = sort(all_corrs);
p_cutoff = sorted_corrs(950);
p_cutoff = 0.45;

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

count_corr = zeros(nROI,nSubj,nCond);
total_corr = zeros(nROI,nSubj,nCond);
for ROI = 1:nROI
	for subj = 1:nSubj
		load([subjList{subj},'Anal/',analName]);
		nVoxels = s.rvf.roi{ROI}.n;
		for condi = 1:nCond
			mean_corr = zeros(1,nVoxels);
			T = 100;
			voxelMeans = zeros(M,nVoxels);
			voxelMeans2 = zeros(M,nVoxels);
			if(condi == 1)
				trainInstances = s.rvf.roi{ROI}.instance.instances(1:M);
			else
				trainInstances = s.rvf.roi{ROI}.instance.instances(M+1:2*M);
			end
			for j = 1:T
				mapping = randperm(27);
				for i = 1:M
					voxelMeans(i,:) = mean(trainInstances{i}(mapping(1:round(N_trials/2)),:));
					voxelMeans2(i,:) = mean(trainInstances{i}(mapping(round(N_trials/2)+1:end),:));
				end
				mean_corr = mean_corr + 1.0/T*diag(corr(voxelMeans,voxelMeans2))';
			end
			count_corr(ROI,subj,condi) = sum(mean_corr > p_cutoff);
			total_corr(ROI,subj,condi) = nVoxels;
		end
	end
end


%%only one
condi = 2;
subj = 3;
ROI = 1;
load([subjList{subj},'Anal/',analName]);
nVoxels = s.lvf.roi{ROI}.n;
mean_corr = zeros(1,nVoxels);
T = 100;
voxelMeans = zeros(M,nVoxels);
voxelMeans2 = zeros(M,nVoxels);
if(condi == 1)
	trainInstances = s.lvf.roi{ROI}.instance.instances(1:M);
else
	trainInstances = s.lvf.roi{ROI}.instance.instances(M+1:2*M);
end
for j = 1:T
	mapping = randperm(27);
	for i = 1:M
		voxelMeans(i,:) = mean(trainInstances{i}(mapping(1:round(N_trials/2)),:));
		voxelMeans2(i,:) = mean(trainInstances{i}(mapping(round(N_trials/2)+1:end),:));
	end
	mean_corr = mean_corr + 1.0/T*diag(corr(voxelMeans,voxelMeans2))';
end
count_corr(ROI,subj,condi) = sum(mean_corr > p_cutoff);
total_corr(ROI,subj,condi) = nVoxels;