% Simulation of the forward encoding model
% --taosheng liu (11/2015)
function simuCinvor(scenario)

switch scenario
    case 1
        disp('Amplitude modulation only');
        kappa=4;
        ampRatio=[1 0.8 0.6 0.4];
        noiseRatio=0.09;
        plotRatio=ampRatio;
        
        m=specifyModel(kappa,'baseModel');
        e=specifyExperiment(m);
        for i=1:length(ampRatio)
            instancesTrain=generateBoldResp(m,e,noiseRatio);
            instancesTest=generateBoldResp(m,e,noiseRatio);
            instancesTrain=cellfun(@(x) mtimes(x,ampRatio(i)),instancesTrain,'UniformOutput',false); %comment this line out to see amplitude effect
            instancesTest=cellfun(@(x) mtimes(x,ampRatio(i)),instancesTest,'UniformOutput',false);
            [channelPrefs thisChannel thisFit thisR2]=fitForwardModel(instancesTrain,e.stimVals,instancesTest,e.stimVals,m,'ampScaleModel');
            allChannelTuning(i,:)=thisChannel;
            allChannelFit(i,:)=thisFit;
            allR2(i)=thisR2;
      end
    case 2
        disp('Noise modulation only');
        kappa=4;
        noiseRatio=[0.09 0.12 0.14 0.18];
        plotRatio=noiseRatio;
        
        m=specifyModel(kappa,'baseModel');
        e=specifyExperiment(m);
        for i=1:length(noiseRatio)
            instancesTrain=generateBoldResp(m,e,noiseRatio(i));
            instancesTest=generateBoldResp(m,e,noiseRatio(i));
            [channelPrefs thisChannel thisFit thisR2]=fitForwardModel(instancesTrain,e.stimVals,instancesTest,e.stimVals,m,'noiseModel');
            allChannelTuning(i,:)=thisChannel;
            allChannelFit(i,:)=thisFit;
            allR2(i)=thisR2;
      end
    case 3
        disp('Tuning width modulation only');
        kappa=[4 3 2 1];
        noiseRatio=[0.09 0.09 0.09 0.09];
        plotRatio=kappa;
        
        for i=1:length(kappa)
            m=specifyModel(kappa(i),'baseModel');
            e=specifyExperiment(m);
            instancesTrain=generateBoldResp(m,e,noiseRatio(i));
            instancesTest=generateBoldResp(m,e,noiseRatio(i));
            [channelPrefs thisChannel thisFit thisR2]=fitForwardModel(instancesTrain,e.stimVals,instancesTest,e.stimVals,m,'kappaModel');
            allChannelTuning(i,:)=thisChannel;
            allChannelFit(i,:)=thisFit;
            allR2(i)=thisR2;
        end
end

h=figure(scenario*10);
set(h,'name',['Scenario ', num2str(scenario)]);
subplot(5,2,1:2);
allChannelFit(:,1)=allChannelFit(:,1)*m.rangeScaleFac;
for i=1:length(plotRatio)
  xp= channelPrefs(1):1:channelPrefs(end)*m.rangeScaleFac;
  yp=fit2VM(allChannelFit(i,1:4),d2r(xp));
  xpOri=xp/m.rangeScaleFac;
  plot(channelPrefs, allChannelTuning(i,:),'o'); hold on
  plot(xpOri,yp,[getcolor(i),'--']);
end
title('Channel tuning functions');

subplot(5,2,3);
bar(r2d(allChannelFit(:,1)));
set(gca,'xticklabel',cellArray(plotRatio));
title('mu');

subplot(5,2,4);
bar(allChannelFit(:,3));
set(gca,'xticklabel',cellArray(plotRatio));
title('baseline');

subplot(5,2,5);
bar(allChannelFit(:,4));
set(gca,'xticklabel',cellArray(plotRatio));
title('Amplitude');

subplot(5,2,6);
bar(allChannelFit(:,5));
set(gca,'xticklabel',cellArray(plotRatio));
title('FWHM');

subplot(5,2,7);
bar(allR2);
set(gca,'xticklabel',cellArray(plotRatio));
title('r2');

subplot(5,2,8);
plot(allR2,allChannelFit(:,3),'o-'); hold on
title('baseline vs. r2')

subplot(5,2,9);
plot(allR2,allChannelFit(:,4),'o-'); hold on
title('amplitude vs. r2');

subplot(5,2,10);
plot(allR2,allChannelFit(:,5),'*--'); 
title('FWHM vs. r2');

return;

% Specify a neuronal model of voxels. Assume neuronal tuning to
% feature ranges from 0 to 360 deg (e.g., direction), and each voxel is made
% up of a mixture of neural tuning functions. In the default setting,
% assume we have 100 voxels, with each voxel containing 180 neurons, which
% are randomly tuning to all possible 360 directions (spaced every 4 deg).
% nvox: number of voxels; neuronsPerVox: number of neurons per voxel
% weighting: distribution of neuronal tuning within a voxel, default to random weighting
% ntuning: number of tuning functions, evenly spaced to cover the 0-360 range
% kappa: concentraion parameter of individual neurons's tuning width, using von Mises function
function m=specifyModel(kappa,figTitle, varargin)
getArgs(varargin, {'nvox=100','neuronsPerVox=180','weighting=random','ntuning=180','plotModel=1'});
totalRangeDeg=180; % Orientation spans 180 deg
rangeScalFac=360/totalRangeDeg;

if rem(neuronsPerVox, ntuning)~=0
  error('We would like to have each voxel contain all possible tuning values');
end
if rem(totalRangeDeg, ntuning)~=0
  error('it would make more sense to have n set of possible tuning functions such as it can be divided by 360');
end

nset=neuronsPerVox/ntuning; %how many sets of the complete neuron tuning is in each voxel
setScaleFac=totalRangeDeg/ntuning; %in case we have fewer possible tuning values
prefStimPerSet= setScaleFac*(0:ntuning-1);
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
m.rangeScaleFac=rangeScalFac;
m.totalRangeDeg=totalRangeDeg;

if plotModel
  smartfig([figTitle,' Left is neural & Right is fMRI'],'reuse');clf
  % show some neuronal tuning curves
  subplot(3,2,1);hold on;
  angle=0:0.01:d2r(totalRangeDeg);
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
  voxIdx=1:nvox/plotVox:nvox;
  plotAngles=0:2:totalRangeDeg; %stimulus values used to generate voxel tuning curve
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
  testAngle=round(rand*totalRangeDeg); %value of the test stimulus, can be any angle
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

% Set up an experiment, with a certain number of stimulus values and number
% of trials per stimulus value, default to 8 and 24, respectively.
function e=specifyExperiment(m,varargin)
getArgs(varargin, {'stimLevel=8','trialPerStim=24'});
e.stimLevel=stimLevel;
e.stimVals=0:m.totalRangeDeg/stimLevel:(m.totalRangeDeg-1);
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
% responses across trials, noise is added. The noiseRatio variable controls the
% magnitude of noise. Also note no time-domain convolution is used as we are not
% modeling fMRI response over time, just the neural responses summed across neurons within a voxel.
function instances=generateBoldResp(m,e,noiseRatio)
instances={};
for i=1:e.stimLevel
  thisStimVal=d2r(e.stimVals(i)); %stimulus value for this current level
  resp=[];
  for j=1:m.neuronsPerVox %loop across neurons within each voxel
    thisPrefStim=d2r(m.prefStim(j));
    resp(:,j)=circ_vmpdf(repmat(thisStimVal,1,m.nvox)*m.rangeScaleFac,thisPrefStim*m.rangeScaleFac,m.kappa);
  end
  wresp=resp.*m.ws; %weight each neurons' response in each voxel
  voxResp=sum(wresp,2); %sum all neurons within a voxel
  voxResp=voxResp'; % reorient to fit the instances structure
  noiseSD=mean(voxResp)*noiseRatio;
  instances{i}=repmat(voxResp, e.trialPerStim, 1)+randn(e.trialPerStim, m.nvox)*noiseSD; %repeat over trials and add noise
end

%fit the forward model to the generated data, first build channels on
%training data, then test channels on test data.
function [channelPrefs, channelTuning, chanFit, r2val]=fitForwardModel(instancesTrain,trainStimVals,instancesTest,testStimVals,m,figTitle, varargin)
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
chanFit=fitTuningWithVM(channelTuning, channelPrefs,m.rangeScaleFac,0);
chanFit=chanFit(1,:);
r2val=r2.overall;


%fit the channel response to a von Mises; note the input chanResp is
%supposed to be two rows, used in real data for ipsi and contra ROI but in
%the context of the simulation they are identical
function fittedVals=fitTuningWithVM(chanResp,channelPref,rangeScaleFac,dispFig)
chanResp=[chanResp;chanResp];

if ~all(chanResp(:,1)==chanResp(:,end))
  disp('First and last element of channel response different, wrapping around it by appending the first element to the end.');
  chanResp=[chanResp,chanResp(:,1)]; %wrap around first and last value
end

paramInit=[pi,1,0,1]; %reasonable guesses of initial value
[paramFit1,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*rangeScaleFac),chanResp(1,:));
[paramFit2,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*rangeScaleFac),chanResp(2,:));

xp= channelPref(1):1:channelPref(end)*rangeScaleFac;
xpOri=xp/rangeScaleFac;

yp1=fit2VM(paramFit1,d2r(xp));
yp1_base0=yp1-min(yp1);
rgHalf1=find(yp1_base0>max(yp1_base0)/2);
fwhm1=xpOri(rgHalf1(end))-xpOri(rgHalf1(1));

yp2=fit2VM(paramFit2,d2r(xp));
yp2_base0=yp2-min(yp2);
rgHalf2=find(yp2_base0>max(yp2_base0)/2);
fwhm2=xpOri(rgHalf2(end))-xpOri(rgHalf2(1));

%convert the mu from 2pi to pi, because it's orientation
paramFit1(1)=paramFit1(1)/rangeScaleFac;
paramFit2(1)=paramFit2(1)/rangeScaleFac;
fittedVals=[[paramFit1, fwhm1];[paramFit2, fwhm2]];

if dispFig %plot the fitted von Mises
  figure;
  subplot(1,2,1);
  plot(channelPref, chanResp(1,:),'o'); hold on
  plot(xpOri,yp1,'r--');
%   yaxis([0 0.5]);
  line([xpOri(rgHalf1(1)) xpOri(rgHalf1(1))], [0 yp1(rgHalf1(1))]);
  line([xpOri(rgHalf1(end)) xpOri(rgHalf1(end))], [0 yp1(rgHalf1(end))]);
  line([xpOri(rgHalf1(1)) xpOri(rgHalf1(end))], [yp1(rgHalf1(1)), yp1(rgHalf1(end))]);
  title(['mu=',num2str(r2d(paramFit1(1))), ' base=',num2str(paramFit1(3)), ' amp=',num2str(paramFit1(4)),' fwhm=',num2str(fwhm1)]);

  subplot(1,2,2);
  plot(channelPref, chanResp(2,:),'o'); hold on
  plot(xpOri,yp2,'r--');
%   yaxis([0 0.5]);
  line([xpOri(rgHalf2(1)) xpOri(rgHalf2(1))], [0 yp2(rgHalf2(1))]);
  line([xpOri(rgHalf2(end)) xpOri(rgHalf2(end))], [0 yp2(rgHalf2(end))]);
  line([xpOri(rgHalf2(1)) xpOri(rgHalf2(end))], [yp2(rgHalf2(1)), yp2(rgHalf2(end))]);
  title(['mu=',num2str(r2d(paramFit2(1))), ' base=',num2str(paramFit2(3)), ' amp=',num2str(paramFit2(4)),' fwhm=',num2str(fwhm2)]);
end



%a 4-parameter circular normal (von Mises) function
%params(1):mean; params(2):kappa; params(3):baseline; params(4):amplitude
function y=fit2VM(params, x)

mu=params(1);
kappa=params(2);
base=params(3);
amp=params(4);

vm=circ_vmpdf(x,mu,kappa);
vm_norm=(vm-min(vm))/(max(vm)-min(vm));
y=base+amp*vm_norm;
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


