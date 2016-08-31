
%simInst.m
%
%
%author : steeve laquitaine
%  date : 160831
%purpose: simulate 2 voxel response instances to be used with the channel
%         encoding model. Both voxel responses increase with the stimulus
%         values input.
%
%usage :
%
%        instances = simInst(1:360/8:360,50)


function instances = simInst(stimValues,numInstances)

figure('color','w')

%number of classes
nClasses = length(stimValues);

%all classes have the same number of instances
numInstances = repmat(numInstances,nClasses,1);

%2 voxel responses clusters
%by class 
%average response cluster by class
mu = 1 : nClasses;

%response dispersion by class 
sigma = eye(nClasses,nClasses);

%structure of response instance x voxel matrices by class
instances = cell(1,nClasses);
for class = 1 : nClasses
    m = mu + (class - 1)*5;
    instances{class} = mvnrnd(m,sigma,numInstances(class))*100;
end

%plot the 2 voxels (x and y axes) 2D response clusters by class (colors)
cols = linspecer(nClasses);
h = nan(nClasses,1);
for i = 1 : nClasses
    hold on; h(i)=plot(instances{i}(:,1),instances{i}(:,2),'.','color',cols(i,:));
    set(h(i),'displayname',num2str(stimValues(i)))
end
legend(h)
xlabel('Voxel 1 responses (AU)')
ylabel('Voxel 2 responses (AU)')



