% getCinvorResponses.m
%
%      usage: instances = getCinvorInstances(m,e)
%         by: justin gardner
%       date: 05/17/16
%    purpose: 
%
function instances = getCinvorInstances(m,e)

% check arguments
if ~any(nargin == [1 2])
  help simulateCinvorResponses
  return
end

% simulate bold responses
% Generate fMRI response, input is the model and experiment (m and e
% structures). For each stimulus value, generate individual trial reponses
% for each voxel. Note the response to a fixed stimulus value would be
% identical across repeats in without any noise. To generate different
% responses across trials, noise is added. The noiseRatio variable controls the
% magnitude of noise. Also note no time-domain convolution is used as we are not
% modeling fMRI response over time, just the neural responses summed across neurons within a voxel.
instances={};
for i=1:e.stimLevel
  thisStimVal=d2r(e.stimVals(i)); %stimulus value for this current level
  resp=[];
  for j=1:m.neuronsPerVox %loop across neurons within each voxel
    thisPrefStim=d2r(m.prefStim(j));
    resp(:,j)=circ_vmpdf(repmat(thisStimVal,1,m.nVoxels)*m.rangeScaleFac,thisPrefStim*m.rangeScaleFac,m.kappa);
  end
  wresp=resp.*m.ws; %weight each neurons' response in each voxel
  voxResp=sum(wresp,2); %sum all neurons within a voxel
  voxResp=voxResp'; % reorient to fit the instances structure
  noiseSD=mean(voxResp)*m.noise;
  instances{i}=repmat(voxResp, e.trialPerStim, 1)+randn(e.trialPerStim, m.nVoxels)*noiseSD; %repeat over trials and add noise
end

% scale by amplitude
instances = cellfun(@(x) mtimes(x,m.amplitude),instances,'UniformOutput',false);


