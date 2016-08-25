% Simulation of the forward encoding model
%
function simuForwardModel

ampIncRatio=2;
ampDecRatio=0.5;

m=specifyModel('baseModel');
e=specifyExperiment;
instancesTrain=generateBoldResp(m,e);
instancesTest=generateBoldResp(m,e);
chParam=fitForwardModel(instancesTrain,e.stimVals,instancesTest,e.stimVals,'baseModel');
keyboard;
% disp('press key to see zscored'); pause;
% fitForwardModel(instancesTrain,e.stimVals,instancesTest,e.stimVals,'zscored','zs=1');

% disp('press key to see amplitude doubled,e.g., train on low contrast, but test on high contrast'); pause;
instancesTest2=cellfun(@(x) mtimes(x,ampIncRatio),instancesTest,'UniformOutput',false);
chParamDouble=fitForwardModel(instancesTrain,e.stimVals,instancesTest2,e.stimVals,['amplitude*',num2str(ampIncRatio)]);

% disp('press key to see amplitude halved,e.g., train on high contrast, but test on low contrast'); pause;
instancesTest3=cellfun(@(x) mtimes(x,ampDecRatio),instancesTest,'UniformOutput',false);
chParamHalf=fitForwardModel(instancesTrain,e.stimVals,instancesTest3,e.stimVals,['amplitude*',num2str(ampDecRatio)]);

% disp('press key to see tuning widened'); pause;
m=specifyModel('widenedTuningModel','kappa=2');
instancesTrain=generateBoldResp(m,e);
instancesTest=generateBoldResp(m,e);
chParamWiden=fitForwardModel(instancesTrain,e.stimVals,instancesTest,e.stimVals,'widenedTuningModel');

h=figure(11);clf
set(h,'name','parameter values');
paramLabels={'mean(rad)','kappa','baseline','amplitude'};
for i=1:4
  subplot(4,1,i);
  thisParam=[chParam(i) chParamDouble(i) chParamHalf(i) chParamWiden(i)];
  bar(thisParam);
  set(gca,'xticklabel',{'original',['Inc',num2str(ampIncRatio)],['Dec',num2str(ampDecRatio)],'widen'});
  title(paramLabels{i});
end
% Specify a neuronal model of voxels. Assume neuronal tuning to
% feature ranges from 0 to 360 deg (e.g., direction), and each voxel is made
% up of a mixture of neural tuning functions. In the default setting,
% assume we have 100 voxels, with each voxel containing 180 neurons, which
% are randomly tuning to all possible 360 directions (spaced every 4 deg).
% nvox: number of voxels; neuronsPerVox: number of neurons per voxel
% weighting: distribution of neuronal tuning within a voxel, default to random weighting
% ntuning: number of tuning functions, evenly spaced to cover the 0-360 range
% kappa: concentraion parameter of individual neurons's tuning width, using von Mises function
function m=specifyModel(figTitle, varargin)
getArgs(varargin, {'nvox=100','neuronsPerVox=180','weighting=random','ntuning=180','kappa=4','plotModel=1'});

if rem(neuronsPerVox, ntuning)~=0
  error('We would like to have each voxel contain all possible tuning values');
end
if rem(360, ntuning)~=0
  error('it would make more sense to have n set of possible tuning functions such as it can be divided by 360');
end

nset=neuronsPerVox/ntuning; %how many sets of the complete neuron tuning is in each voxel
scaleFac=360/ntuning; %in case we have fewer possible tuning values
prefStimPerSet= scaleFac*(0:ntuning-1);
prefStim=[];
for stim=1:nset
  prefStim=[prefStim,prefStimPerSet];
end

%generate weights for each voxel, returns a nvox by neuronsPerVox matrix
if strcmp(weighting,'random')
  for i=1:nvox
    thisW=rand(1,neuronsPerVox); %for each voxel, generate random weights for each tuning function
    thisW=thisW/sum(thisW); %normalize these weights, so that they sum to 1
    ws(i,:)=thisW;
  end
else
  error('other weighting scheme is not implemented yet');
end

%put these parameters in the m structure
m.ws=ws;
m.weighting=weighting;
m.nvox=nvox;
m.neuronsPerVox=neuronsPerVox;
m.ntuning=ntuning;
m.prefStim=prefStim;
m.kappa=kappa;

if plotModel
  smartfig([figTitle,' Left is neural & Right is fMRI'],'reuse');clf
  % show some neuronal tuning curves
  subplot(3,2,1);hold on;
  angle=0:0.01:2*pi;
  int=neuronsPerVox/10; %select every 10 neurons
  for i=1:int:neuronsPerVox
    thisNeuralTuning=circ_vmpdf(angle,d2r(prefStim(i)), kappa);
    plot(r2d(angle),thisNeuralTuning,getcolor(ceil(i/int)));
  end
  ylabel('Response');
  xlabel('Stim values');
  title('Tuning curves of 10 representative neurons');
  
  % show some voxel tuning curves
  subplot(3,2,3); hold on
  plotVox=10; %select 10 voxels to plot
  voxIdx=1:nvox/plotVox:nvox;
  plotAngles=0:2:358; %stimulus values used to generate voxel tuning curve
  for istim=1:length(plotAngles)
    thisStimVal=d2r(plotAngles(istim));
    resp=[];
    %for each voxel, get all neuronal responses for the current stimulus
    for j=1:neuronsPerVox
      thisPrefStim=d2r(prefStim(j));
      resp(:,j)=circ_vmpdf(repmat(thisStimVal,1,plotVox), thisPrefStim, kappa);
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
  testAngle=round(rand*360); %value of the test stimulus, can be any angle
  %for each neuronal tuning channel, caculate its response to the test stimulus
  for i=1:length(prefStimPerSet)
    popResp(i)=circ_vmpdf(d2r(testAngle), d2r(prefStimPerSet(i)), kappa);
  end
  plot(prefStimPerSet, popResp, 'o-','lineWidth',2);
  ylabel('Response');
  xlabel('Preferred stim');
  vline(testAngle);
  title(['Neuronal population response to a stimulus at ', num2str(testAngle)]);
  
end

% Set up an experiment, with a certain number of stimulus values and number
% of trials per stimulus value, default to 8 and 24, respectively.
function e=specifyExperiment(varargin)
getArgs(varargin, {'stimLevel=8','trialPerStim=24'});
e.stimLevel=stimLevel;
e.stimVals=0:360/stimLevel:359;
e.trialPerStim=trialPerStim;
allStimVal=[];
for i=1:stimLevel
  allStimVal=[allStimVal, repmat(e.stimVals(i),1,trialPerStim)];
end
e.allStimVal=allStimVal; %contains stimulus values on each trial

% Generate fMRI response, input is the model and experiment (m and e
% structures). For each stimulus value, generate individual trial reponses
% for each voxel. Note the response to a fixed stimulus value would be
% identical across repeats in without any noise. To generate different
% responses across trials, noise is added. The sdratio variable controls the
% magnitude of noise. Also note no signal convolution is used as we are not
% modeling fMRI response over time, just the neural responses summed across neurons within a voxel.
function instances=generateBoldResp(m,e,varargin)
getArgs(varargin, {'sdratio=4'}); %assume noise SD is a quarter of mean signal
instances={};
for i=1:e.stimLevel
  thisStimVal=d2r(e.stimVals(i)); %stimulus value for this current level
  resp=[];
  for j=1:m.neuronsPerVox %loop across neurons within each voxel
    thisPrefStim=d2r(m.prefStim(j));
    resp(:,j)=circ_vmpdf(repmat(thisStimVal,1,m.nvox),thisPrefStim,m.kappa);
  end
  wresp=resp.*m.ws; %weight each neurons' response in each voxel
  voxResp=sum(wresp,2); %sum all neurons within a voxel
  voxResp=voxResp'; % reorient to fit the instances structure
  noiseSD=mean(voxResp)/sdratio;
  instances{i}=repmat(voxResp, e.trialPerStim, 1)+randn(e.trialPerStim, m.nvox)*noiseSD; %repeat over trials and add noise
end

%fit the forward model to the generated data, first build channels on
%training data, then test channels on test data.
function paramFit=fitForwardModel(instancesTrain,trainStimVals,instancesTest,testStimVals,figTitle, varargin)
getArgs(varargin,{'zs=0'});

plotMeanInstances(instancesTrain, instancesTest,trainStimVals,testStimVals);

if zs
  channel=buildChannels(instancesTrain,trainStimVals,'dispChannels=1','zscore=1');
else
  channel=buildChannels(instancesTrain,trainStimVals,'dispChannels=1');
end

[channelTuning, r2, classifyCorrTotal, stimValVector, predStimVal]=testChannels(instancesTest,testStimVals,channel);
channelPrefs=[channel.channelPref, channel.channelPref(end)+unique(diff(channel.channelPref))]; %just wrap around the end
channelTuning=[channelTuning, channelTuning(1)]; %repeat the first response to wrap around
%figure out the ideal response, basically the ideal response to the center stimulus
centerIdx=find(channel.channelPref==channel.span/2);
idealResp=channel.idealStimResponse(centerIdx,:);
idealResp=[idealResp, idealResp(1)];

smartfig([figTitle,' Left is neural & Right is fMRI'],'reuse');
subplot(3,2,2);
plot(channelPrefs,channelTuning,'o-','linewidth',2); hold on
plot(channelPrefs,idealResp,'g');
title('Channel tuning');
legend({'Avg test resp','Ideal resp'});

subplot(3,2,4);
bar([r2.clAcc;r2.clOneshot;r2.clAvg]');
yaxis([0 1]);
title(['r2 for each stim class, grand avg r2:',num2str(r2.overall)]);
set(gca,'xticklabel',num2cell(trainStimVals));
legend({'accu','oneshot','avg'});

subplot(3,2,6);
% predErr=r2d(circ_dist(d2r(stimValVector),d2r(predStimVal)));
% hist(predErr);
% title(['Prediction errors, accuracy=',num2str(classifyCorrTotal(1)/classifyCorrTotal(2))]);
plot([r2.voxAcc;r2.voxOneshot;r2.voxAvg]');
legend({'accu','oneshot','avg'});
title('r2 value for each voxel');
xlabel('voxel number');
% fit the channel response to a von Mises
paramInit=[pi,1,0,1]; %reasonable guesses of initial value
[paramFit,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPrefs),channelTuning);
% disp(['Fitted params of channel tuning, mean=',num2str(r2d(paramFit(1))),' kappa=',num2str(paramFit(2)),...
%   ' baseline=',num2str(paramFit(3)),' amplitude=',num2str(paramFit(4))]);
%plot the fitted von Mises
xp=channelPrefs(1):1:channelPrefs(end);
yp=fit2VM(paramFit,d2r(xp));
subplot(3,2,2);
plot(xp,yp,'--');
legend({'Avg test resp','Ideal resp','Fitted resp'});


% %get r2 values for predicted vs. actual test instances
% function thisClassR2=getR2PredInstances(predInstances, instancesTest)
% for i=1:length(instancesTest)
%   thisTestClass=instancesTest{i};
%   for j=1:size(thisTestClass,1)
%     thisTestTrial=thisTestClass(j,:);
%     sumOfSquaresResidual = sum((thisTestTrial-predInstances{i}).^2);
%     thisR2(j)=1-sumOfSquaresResidual./sum(thisTestTrial.^2);
%   end
%   thisClassR2(i)=mean(thisR2);
% end

%a 4-parameter circular normal (von Mises) function
function y=fit2VM(params, x)

mu=params(1);
kappa=params(2);
base=params(3);
amp=params(4);

y=base+amp*circ_vmpdf(x,mu,kappa);
y=y'; %circ_vmpdf always return column vector, turn it into row

%Just plot the average instances for each stim class
function plotMeanInstances(instancesTrain, instancesTest, trainStimVals,testStimVals)

temp=cellfun(@(x) mean(x(:)),instancesTrain,'UniformOutput',false);
meanInsTrain=cell2mat(temp);

temp=cellfun(@(x) mean(x(:)),instancesTest,'UniformOutput',false);
meanInsTest=cell2mat(temp);

smartfig('Mean Instances of each class','reuse');clf
hold on
plot(trainStimVals, meanInsTrain,'go-');
plot(testStimVals, meanInsTest,'rs-');
legend({'train','test'});
yaxis([0,0.5]);


