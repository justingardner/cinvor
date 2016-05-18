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
%
function m = setCinvorModel(e,varargin)

% get arguments
getArgs(varargin, {'kappa=4','amplitude=1','noise=0.1','nVoxels=100','neuronsPerVox=180','weighting=random','nTuning=180'});

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

%generate weights for each voxel, returns a nVoxels by neuronsPerVox matrix
if strcmp(weighting,'random')
  for i=1:nVoxels
    thisW=rand(1,neuronsPerVox); %for each voxel, generate random weights for each tuning function
    thisW=thisW/sum(thisW); %normalize these weights, so that they sum to 1
    ws(i,:)=thisW;
  end
else
  disp('(setCinvorModel) Other weighting scheme is not implemented yet');
end

%put these parameters in the m structure
m.ws=ws;
m.weighting=weighting;
m.nVoxels=nVoxels;
m.neuronsPerVox=neuronsPerVox;
m.nTuning=nTuning;
m.prefStim=prefStim;
m.kappa=kappa;
m.amplitude = amplitude;
m.noise = noise;
m.rangeScaleFac=rangeScalFac;

if e.plotModel
  smartfig([e.figTitle,' Left is neural & Right is fMRI'],'reuse');clf
  % show some neuronal tuning curves
  subplot(3,2,1);hold on;
  angle=0:0.01:d2r(e.totalRangeDeg);
  int=neuronsPerVox/9; %select 9 neurons to plot
  for i=1:int:neuronsPerVox
    thisNeuralTuning=circ_vmpdf(angle*rangeScalFac,d2r(prefStim(i))*rangeScalFac, kappa);
    plot(r2d(angle),thisNeuralTuning,getcolor(ceil(i/int)));
  end
  ylabel('Response');
  xlabel('Stim values');
  title('Tuning curves of 10 representative neurons');
  
  % show some voxel tuning curves
  subplot(3,2,3); hold on
  plotVox=10; %select 10 voxels to plot
  voxIdx=1:nVoxels/plotVox:nVoxels;
  plotAngles=0:2:e.totalRangeDeg; %stimulus values used to generate voxel tuning curve
  for istim=1:length(plotAngles)
    thisStimVal=d2r(plotAngles(istim));
    resp=[];
    %for each voxel, get all neuronal responses for the current stimulus
    for j=1:neuronsPerVox
      thisPrefStim=d2r(prefStim(j));
      resp(:,j)=circ_vmpdf(repmat(thisStimVal,1,plotVox)*rangeScalFac, thisPrefStim*rangeScalFac, kappa);
    end
    wresp=resp.*ws(voxIdx,:); %weight the neuronal response with pre-generated weights
    voxResp(:,istim)=sum(wresp,2); %sum across all neurons within a voxel
  end
  for v=1:plotVox
    plot(plotAngles, voxResp(v,:), getcolor(v));
  end
  ylabel('Response');
  xlabel('Stim values');
  title('Tuning profile for 10 randomly selected voxels');
  
  %for a single stimulus, plot the neuronal population response
  subplot(3,2,5); hold on
  testAngle=round(rand*e.totalRangeDeg); %value of the test stimulus, can be any angle
  %for each neuronal tuning channel, caculate its response to the test stimulus
  for i=1:length(prefStimPerSet)
    popResp(i)=circ_vmpdf(d2r(testAngle)*rangeScalFac, d2r(prefStimPerSet(i))*rangeScalFac, kappa);
  end
  plot(prefStimPerSet, popResp, 'o-','lineWidth',2);
  ylabel('Response');
  xlabel('Preferred stim');
  vline(testAngle);
  title(['Neuronal population response to a stimulus at ', num2str(testAngle)]);
  
end

