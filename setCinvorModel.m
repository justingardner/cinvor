% setCinvorModel.m
%
%      usage: setCinvorModel()
%         by: Taosheng Liu
%    purpose: 
%
% Specify a neuronal model of voxels. Assume neuronal tuning to
% feature ranges from 0 to 360 deg (e.g., direction), and each voxel is made
% up of a mixture of neural tuning functions. In the default setting,
% assume we have 100 voxels, with each voxel containing 180 neurons, which
% are randomly tuning to all possible 360 directions (spaced every 4 deg).
% nvox: number of voxels; neuronsPerVox: number of neurons per voxel
% weighting: distribution of neuronal tuning within a voxel, default to random weighting
% ntuning: number of tuning functions, evenly spaced to cover the 0-360 range
% kappa: concentraion parameter of individual neurons's tuning width, using von Mises function
% myWeights: fixed selection of weights (must use weighting=fixed)
% alpha: proportion of voxels containing information
function m = setCinvorModel(e,varargin)

% get arguments %%%%%%kappa=2.5, noise = 8
getArgs(varargin, {'kappa=4','amplitude=1','noise=10','nVoxels=78','neuronsPerVox=180','weighting=random','nTuning=180','myWeights=0','uniformFactor=0','tao=.1','exponent=2','alphaParam=1','kConcentrated=[]'});

% range of orientations span 180
rangeScalFac=360/e.totalRangeDeg;

% sanity checks
if rem(neuronsPerVox, nTuning)~=0
  disp('(setCinvorModel) We would like to have each voxel contain all possible tuning values');
  keyboard
end
if rem(e.totalRangeDeg, nTuning)~=0
  disp('(setCinvorModel) It would make more sense to have n set of possible tuning functions such as it can be divided by 360');
  keyboard
end

% how many sets of the complete neuron tuning is in each voxel
nset=neuronsPerVox/nTuning;
%in case we have fewer possible tuning values
setScaleFac=e.totalRangeDeg/nTuning;
prefStimPerSet= setScaleFac*(0:nTuning-1);
prefStim=[];
for stim=1:nset
  prefStim=[prefStim,prefStimPerSet];
end

% angles over which tuning and weights are computed
angles = 0:e.totalRangeDeg-1;

% initialze weights
neuronVoxelWeights = zeros(nVoxels,neuronsPerVox);
%generate weights for each voxel, returns a nVoxels by neuronsPerVox matrix
if strcmp(weighting,'random')
  % cycle over voxels
  for iVoxel = 1:nVoxels
    if(rand() < alphaParam)
      % get uniform random weighting of neurons - if uniformFactor is set, then adds
      % weighting of all neurons to every voxel
      neuronVoxelWeights(iVoxel,:) = uniformFactor*ones(1,neuronsPerVox)+rand(1,neuronsPerVox); %for each voxel, generate random weights for each tuning function
    else
      neuronVoxelWeights(iVoxel,:) = zeros(1,neuronsPerVox);
    end
  end
elseif strcmp(weighting,'fixed')
  % passed in weight matrix
  weights = myWeights;
elseif strcmp(weighting,'concentrated')
  % this picks weights for each voxel conentrated around
  % a random orientation for each voxel - according to the
  % von mises distribution with concentration parameter kConcentrated
  % first pick a random orientation for each voxel
  m.voxelPreferredOrientation = floor(rand(1,nVoxels)*180);
  if isinf(kConcentrated) & isinf(kappa)
    % make sure the stimulus values are in the set, otherwise the simulation will
    % have some stimulus creating no response.
    while ~isempty(setdiff(round(e.stimVals),m.voxelPreferredOrientation))
      m.voxelPreferredOrientation = floor(rand(1,nVoxels)*180);
    end
  end
  % get a canonical von mises function which we circular shift below
  % to get each preferred orientation (this is done for speed).
  vonMises = fitVonMises(angles,[],0,kConcentrated,180);
  for iVoxel = 1:nVoxels
    if isinf(kConcentrated)
      % if concentration parameter is inf then it means to use a stick function
      neuronVoxelWeights(iVoxel,1:length(angles)) = 0;
      neuronVoxelWeights(iVoxel,mod(m.voxelPreferredOrientation(iVoxel),180)+1) = 1;
    else
      % get von mises distribution centered around preferred orientation for voxel (Wrapping at 180 instead of 360)
      neuronVoxelWeights(iVoxel,:) = circshift(vonMises,mod(m.voxelPreferredOrientation(iVoxel),180));
    end
    % and normalize weights to 1 (not sure this is needed, since all voxels will have the same weights)
    neuronVoxelWeights(iVoxel,:) = neuronVoxelWeights(iVoxel,:)./sum(neuronVoxelWeights(iVoxel,:));
  end
else
  disp('(setCinvorModel) Other weighting scheme is not implemented yet');
end

% get a von mises function (which is used below for every possible tuning, by 
% circular shifting
vonMises = fitVonMises(angles,[],0,kappa,180);

% compute tuning of neurons
for iNeuron = 1:nTuning
  % if concentration parameter is inf then it means to use a stick function
  if isinf(kappa)
    m.neuralTuning(iNeuron,1:length(angles)) = 0;
    m.neuralTuning(iNeuron,find(round(prefStim(iNeuron))==angles)) = 1;
  else
    % compute neural response as a vonmises function (wrapping at 180 instead of 360)
    m.neuralTuning(iNeuron,:) = circshift(vonMises,iNeuron-1);
    m.neuralTuning(iNeuron,:) = m.neuralTuning(iNeuron,:)./max(m.neuralTuning(iNeuron,:));
  end
end

% cycle over stimuli and calculate response of neurons
for iStimVal = 1:e.stimLevel
  % initialize neural response matrix
  neuralResponse = zeros(nVoxels,neuronsPerVox);
  % loop across neurons, computing their response to the stimulus
  for jNeuron=1:neuronsPerVox 
    if isinf(kappa)
      % if we have stick functions then respond if there is a match
      % note that we round here to deal with stimulus values like 22.5 which
      % would otherwise (given the degree spacing of tuning never elicit a response.
      if round(e.stimVals(iStimVal)) == round(prefStim(jNeuron))
	neuralResponse(:,jNeuron) = 1;
      else
	neuralResponse(:,jNeuron) = 0;
      end
    else
      % preferred stimulus for each neuron
      neuralResponse(:,jNeuron)= fitVonMises(e.stimVals(iStimVal),[],prefStim(jNeuron),kappa,180);
    end
  end
  % weight each neurons' response in each voxel
  weightedResponse = neuralResponse .* neuronVoxelWeights; 
  %sum all neurons within a voxel
  voxelResponse = sum(weightedResponse,2)'; 
  % put into instances
  m.voxelResponse{iStimVal}=repmat(voxelResponse*amplitude, e.trialPerStim, 1);
end

%put these parameters in the m structure
m.neuronVoxelWeights = neuronVoxelWeights;
m.weighting=weighting;
m.nVoxels=nVoxels;
m.neuronsPerVox=neuronsPerVox;
m.nTuning=nTuning;
m.prefStim=prefStim;
m.kappa=kappa;
m.amplitude = amplitude;
m.noiseSTD = noise;
m.rangeScaleFac=rangeScalFac;
m.uniformFactor = uniformFactor;
m.tao = tao;
m.exponent = exponent;
m.alphaParam = alphaParam;
m.kConcentrated = kConcentrated;

% testing plots
doTestPlot = 0;
if doTestPlot
  testPlot(m,e,neuralResponse)
end


%%%%%%%%%%%%%%%%%%
%    testPlot    %
%%%%%%%%%%%%%%%%%%
function testPlot(m,e,neuralResponse)

mlrSmartfig('setCinvorModel','reuse');clf;
subplot(2,3,4);
plot(m.neuronVoxelWeights');
title(sprintf('All %i voxel weights',m.nVoxels));
subplot(2,3,1);
randVox = ceil(rand*m.nVoxels);
plot(m.neuronVoxelWeights(randVox,:));
hold on
vline(m.voxelPreferredOrientation(randVox)+1);
vline(mod(m.voxelPreferredOrientation(randVox)+fitVonMises([],[],'kappa',m.kConcentrated),180)+1);;
vline(mod(m.voxelPreferredOrientation(randVox)-fitVonMises([],[],'kappa',m.kConcentrated),180)+1);;
title(sprintf('Vox: %i pref: %i halfWidth: %0.1f',randVox,m.voxelPreferredOrientation(randVox),fitVonMises([],[],'kappa',m.kConcentrated)));
hline(max(m.neuronVoxelWeights(randVox,:))/2);
% plot neural tuning
subplot(2,3,5);
plot(m.neuralTuning');
title(sprintf('All %i neurons tuning',m.neuronsPerVox));
subplot(2,3,2);
randNeuron = ceil(rand*m.neuronsPerVox);
plot(m.neuralTuning(randNeuron,:));
hold on
vline(randNeuron);
vline(mod(randNeuron+fitVonMises([],[],'kappa',m.kappa),180));;
vline(mod(randNeuron-fitVonMises([],[],'kappa',m.kappa),180));;
title(sprintf('Neuron: %i pref: %i halfWidth: %0.1f',randNeuron,randNeuron,fitVonMises([],[],'kappa',m.kappa)));
hline(max(m.neuralTuning(randVox,:))/2);
% plot last stimulus response
subplot(2,3,3);
plot(neuralResponse');
title(sprintf('Neural response'));
hold on
vline(e.stimVals(end)+1);
subplot(2,3,6);
plot(m.voxelResponse{end}');
title(sprintf('Voxel response'));
