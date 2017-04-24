% fitCinvorForwardModel.m
%
%      usage: fitCinvorForwardModel()
%         by: Taosheng Liu
%       date: 11/2015
%    purpose: fit the forward model to the generated data, first build channels on
%             training data, then test channels on test data.
%
%
function [channelPrefs, channelTuning, chanFit, r2val] = fitCinvorForwardModel(instancesTrain,trainStimVals,instancesTest,testStimVals,m,e,varargin)
dispFig = 0;
% get arguments
getArgs(varargin,{'zs=0','model=sinFilter','exponent=7'});
if(dispFig)
  plotMeanInstances(instancesTrain, instancesTest,trainStimVals,testStimVals);
end

if zs
  channel=buildChannels(instancesTrain,trainStimVals,'dispChannels=0','zscore=1','fitNoise',0,'model',model,'exponent',exponent);
else
  channel=buildChannels(instancesTrain,trainStimVals,'dispChannels=0','fitNoise',0,'model',model,'exponent',exponent);
end

[channelTuning, r2, classifyCorrTotal, stimValVector, predStimVal]=testChannels(instancesTest,testStimVals,channel,'fitNoise',0);
channelPrefs=[channel.channelPref, channel.channelPref(end)+unique(diff(channel.channelPref))]; %just wrap around the end
channelTuning=[channelTuning, channelTuning(1)]; %repeat the first response to wrap around
%figure out the ideal response, basically the ideal response to the center stimulus
centerIdx=find(channel.channelPref==channel.span/2);
idealResp=channel.idealStimResponse(centerIdx,:);
idealResp=[idealResp, idealResp(1)];
if(dispFig)
  mlrSmartfig(e.figTitle,'reuse'); 
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
end
chanFit=fitTuningWithVM(channelTuning, channelPrefs,m.rangeScaleFac,0);
chanFit=chanFit(1,:);
r2val=mean(r2.voxOneshot);

%fit the channel response to a von Mises; note the input chanResp is
%supposed to be two rows, used in real data for ipsi and contra ROI but in
%the context of the simulation they are identical
function fittedVals=fitTuningWithVM(chanResp,channelPref,rangeScaleFac,dispFig)
chanResp=[chanResp;chanResp];
fixedMean =1;

if ~all(chanResp(:,1)==chanResp(:,end))
  disp('First and last element of channel response different, wrapping around it by appending the first element to the end.');
  chanResp=[chanResp,chanResp(:,1)]; %wrap around first and last value
end

paramInit=[pi,1,0,1]; %reasonable guesses of initial value
if(fixedMean)
  [paramFit1,resnorm,residual,exitflag]=lsqcurvefit(@vonMises,paramInit,d2r(channelPref*rangeScaleFac),chanResp(1,:),[pi-0.01,0,0,0],[pi+0.01,1000,1000,1000]);
  [paramFit2,resnorm,residual,exitflag]=lsqcurvefit(@vonMises,paramInit,d2r(channelPref*rangeScaleFac),chanResp(2,:),[pi-0.01,0,0,0],[pi+0.01,1000,1000,1000]);
else
  [paramFit1,resnorm,residual,exitflag]=lsqcurvefit(@vonMises,paramInit,d2r(channelPref*rangeScaleFac),chanResp(1,:));
  [paramFit2,resnorm,residual,exitflag]=lsqcurvefit(@vonMises,paramInit,d2r(channelPref*rangeScaleFac),chanResp(2,:));
end
xp= channelPref(1):1:channelPref(end)*rangeScaleFac;
xpOri=xp/rangeScaleFac;

yp1=vonMises(paramFit1,d2r(xp));
yp1_base0=yp1-min(yp1);
rgHalf1=find(yp1_base0>max(yp1_base0)/2);
fwhm1=xpOri(rgHalf1(end))-xpOri(rgHalf1(1));

yp2=vonMises(paramFit2,d2r(xp));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    plotMeanInstances    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMeanInstances(instancesTrain, instancesTest, trainStimVals,testStimVals)

%Just plot the average instances for each stim class
temp=cellfun(@(x) mean(x(:)),instancesTrain,'UniformOutput',false);
meanInsTrain=cell2mat(temp);

temp=cellfun(@(x) mean(x(:)),instancesTest,'UniformOutput',false);
meanInsTest=cell2mat(temp);

mlrSmartfig('Mean Instances of each class','reuse');clf
hold on
plot(trainStimVals, meanInsTrain,'go-');
plot(testStimVals, meanInsTest,'rs-');
legend({'train','test'});
xlabel('Orientations');
ylabel('response mean');
%yaxis([0,0.5]);

