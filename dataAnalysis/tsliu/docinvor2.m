%% docinvor.m
%
%        $Id:$
%      usage: docinvor2(rebuild)
%         by: taosheng liu
%       date: 07/25/14
%    purpose: Analyze cinvor2 experiment
%    rebuild=1 -- redo the analysis; rebuild=0 -- use exisiting analysis;
%    rebuild=2 -- extract all instances and save it to disk (for Justin)
%%
function docinvor2(rebuild)
% % BB: build on both contrast and test on low and high separately
%   docinvorOnce('meanInsBB',['rebuild=',num2str(rebuild)]); 
%   docinvorOnce('deconInsBB','instanceMethod=deconv',['rebuild=',num2str(rebuild)],'canonicalType=all'); 
%   docinvorOnce('decon1gInsBB','instanceMethod=deconv',['rebuild=',num2str(rebuild)]);
% % disp('##########################PRESS ANY KEY TO CONTINUE############################'); pause;

% % Standard: build on high contrast and test on high, build on low and test on low, separately
%   docinvorOnce('meanIns',['rebuild=',num2str(rebuild)]); return;
%   docinvorOnce('deconIns','instanceMethod=deconv',['rebuild=',num2str(rebuild)],'canonicalType=all');
docinvorOnce('decon1gIns','instanceMethod=deconv',['rebuild=',num2str(rebuild)]);
 
% % XC: cross contrast, build on low and test on high, and vice versa
%   docinvorOnce('meanInsXC',['rebuild=',num2str(rebuild)]); %use mean instance method
%   docinvorOnce('deconInsXC','instanceMethod=deconv',['rebuild=',num2str(rebuild)],'canonicalType=all'); %decon instance, do not fit with smooth gamma
%   docinvorOnce('decon1gInsXC','instanceMethod=deconv',['rebuild=',num2str(rebuild)]); %decon instance, fit with smooth 1gamma



%% do classification and forward modeling with one setting
function docinvorOnce(saveName, varargin)

getArgs(varargin,{'instanceMethod=mean','rebuild=1','dispFig=0','verbose=1','doClassification=1','doForwardModel=1','canonicalType=[]'});
mrQuit(0);
% pack up into status structure
s.instanceMethod = instanceMethod;
s.doClassification = doClassification;
s.doForwardModel = doForwardModel;
s.rebuild = rebuild;
s.dispFig = dispFig;
s.verbose = verbose;
s.canonicalType=canonicalType;
s.args = varargin;
s.saveName=saveName;

if ~isempty(strfind(saveName,'XC'))
  disp('Do cross-contrast validation, build on high test on low and vice versa');
  s.whichCond=0;
  s = doCrossContrastValidation(s);
elseif ~isempty(strfind(saveName,'BB'))
  disp('Do cross-contrast build on both and test on both');
  s.whichCond=3;
  s.mode='BB';
  s = doLeaveOneOut(s);
else
  disp('Do within-contrast cross-run validation');
  s.whichCond=3; %note 3 does orientation decoding/forward modeling for low and high contrast separately, can change this value to 1 or 2 see doLeaveOneOut for more info
  s.mode='WI';
  s = doLeaveOneOut(s); %this runs classification and forward modeling within contrast
end

%% train on low test on high and vice versa
function s = doCrossContrastValidation(s)
roiList={'lV1res','rV1res','lV1','rV1'};
if s.rebuild
  lvfStr='contrast1 _x_ orientation1';
  rvfStr='contrast2 _x_ orientation2';
  s.lvf=crossContrastValidation(s,roiList,lvfStr,'orientation1','canonicalType',s.canonicalType); 
  s.rvf=crossContrastValidation(s,roiList,rvfStr,'orientation2','canonicalType',s.canonicalType);
  s.roiList=roiList;
  disp(['Saving result: ', 'Anal/',s.saveName]);
  save(['Anal/',s.saveName], 's');
else
  disp(['Loading previously saved analysis: ', pwd, '/', s.saveName]);
  load(['Anal/',s.saveName]);
  roiList=s.roiList;
end

fields=fieldnames(s.lvf);
fmIdx=find(strncmp(fields,'fm',2));
for i=1:length(fmIdx) %TSL: note this combining of left and right VF result is a hack, should improve in the future (match ROI names etc).
  thisFM=fields{fmIdx(i)};
  s.bvf.(thisFM).chanResp=(s.lvf.(thisFM).chanResp + s.rvf.(thisFM).chanResp([2,1],:))/2;
  s.bvf.(thisFM).chanClass=(s.lvf.(thisFM).chanClass + s.rvf.(thisFM).chanClass([2,1],:))/2;
  s.bvf.(thisFM).meanIns=(s.lvf.(thisFM).meanIns + s.rvf.(thisFM).meanIns([2,1],:))/2;
end
plotChannelResultXC(s,roiList,1/8,s.saveName);
disp(['Re-Saving result: ', 'Anal/',s.saveName]);
save(['Anal/',s.saveName], 's','-v7.3');

%% leave one run out, either within contrast or train on both and test on both
function s=doLeaveOneOut(s) 
roiList={'lV1res','rV1res','lV1','rV1'}; 

switch s.whichCond
  case 1
    disp('%%%%%%%%%%%%%%%%%%Contrast only%%%%%%%%%%%%%%%%%%');
    lvfStr='contrast1';
    rvfStr='contrast2';
  case 2
    disp('%%%%%%%%%%%%%%%%%%Orientation only%%%%%%%%%%%%%%%%%%');
    lvfStr='orientation1';
    rvfStr='orientation2';
  case 3
    disp('%%%%%%%%%%%%%%%%%%Both Contrast and Orientation %%%%%%%%%%%%%%%%%%');
    lvfStr='contrast1 _x_ orientation1';
    rvfStr='contrast2 _x_ orientation2';
end

if s.rebuild
  if strcmp(s.mode,'WI')
    s.lvf=leaveOneRunOut(s,roiList,lvfStr,'orientation1','canonicalType',s.canonicalType);
    s.rvf=leaveOneRunOut(s,roiList,rvfStr,'orientation2','canonicalType',s.canonicalType);
  elseif strcmp(s.mode,'BB')
    s.lvf=leaveOneRunOut(s,roiList,lvfStr,'orientation1','canonicalType',s.canonicalType);
    s.rvf=leaveOneRunOut(s,roiList,rvfStr,'orientation2','canonicalType',s.canonicalType);
  end
  if s.rebuild==1
    disp(['Saving result: ', 'Anal/',s.saveName]);
    save(['Anal/',s.saveName], 's');
  elseif s.rebuild==2
    lvf=s.lvf;
    rvf=s.rvf;
    previousName=s.saveName;
    clear s;
    disp(['Loading previously saved analysis: ', previousName]);
    load(['Anal/',previousName]);
    s.lvf.roi=lvf.roi;
    s.rvf.roi=rvf.roi;
    disp(['Re-saving result with ROI instances added: ', 'Anal/',s.saveName]);
    save(['Anal/',s.saveName], 's');
    return;
  end

else
  disp(['Loading previously saved analysis: ', s.saveName]);
  load(['Anal/',s.saveName]);
  roiList={s.lvf.roi{1}.name,s.lvf.roi{2}.name};
end

chance=1/8; %1/length(s.lvf.stimvol); %only used to plot a hline at chance
if s.doClassification
  fields=fieldnames(s.lvf.pcorrLow);
  for i=1:length(fields)
    s.bvf.pcorrLow.(fields{i})=(s.lvf.pcorrLow.(fields{i}) + fliplr(s.rvf.pcorrLow.(fields{i})))/2;
    s.bvf.pcorrHigh.(fields{i})=(s.lvf.pcorrHigh.(fields{i}) + fliplr(s.rvf.pcorrHigh.(fields{i})))/2;
  end
  plotAllAccu(s,roiList,chance,s.saveName,'Low');
  plotAllAccu(s,roiList,chance,s.saveName,'High');
end

if s.doForwardModel
  fields=fieldnames(s.lvf.fmLow.chanResp);
%   fmIdx=find(strncmp(fields,'fm',2));
  for i=1:length(fields)
%     thisFM=fields{fmIdx(i)};
    s.bvf.fmLow.chanResp.(fields{i})=(s.lvf.fmLow.chanResp.(fields{i}) + flipThis(s.rvf.fmLow.chanResp.(fields{i})))/2;
    s.bvf.fmLow.chanClass.(fields{i})=(s.lvf.fmLow.chanClass.(fields{i}) + flipThis(s.rvf.fmLow.chanClass.(fields{i})))/2;
    s.bvf.fmHigh.chanResp.(fields{i})=(s.lvf.fmHigh.chanResp.(fields{i}) + flipThis(s.rvf.fmHigh.chanResp.(fields{i})))/2;
    s.bvf.fmHigh.chanClass.(fields{i})=(s.lvf.fmHigh.chanClass.(fields{i}) + flipThis(s.rvf.fmHigh.chanClass.(fields{i})))/2;
    s.bvf.fmLow.avgR2.(fields{i})=averageIpsiContra(s.lvf.fmLow.avgR2.(fields{i}), s.rvf.fmLow.avgR2.(fields{i}));
    s.bvf.fmHigh.avgR2.(fields{i})=averageIpsiContra(s.lvf.fmHigh.avgR2.(fields{i}), s.rvf.fmHigh.avgR2.(fields{i}));
    s.bvf.fmLow.meanIns.(fields{i})=(s.lvf.fmLow.meanIns.(fields{i}) + s.rvf.fmLow.meanIns.(fields{i})([2,1],:))/2;
    s.bvf.fmHigh.meanIns.(fields{i})=(s.lvf.fmHigh.meanIns.(fields{i}) + s.rvf.fmHigh.meanIns.(fields{i})([2,1],:))/2;
    plotChannelR2(s,fields{i},roiList,s.saveName);
    plotChannelResult(s,fields{i},roiList,chance,s.saveName);
  end
end
disp(['Re-Saving result: ', 'Anal/',s.saveName]);
save(['Anal/',s.saveName], 's');

%% average values for ipis and contra hemifield
function bvf=averageIpsiContra(lvf, rvf)
fields=fieldnames(lvf);
for i=1:length(fields)
  if ~isempty(strfind(fields{i},'vox'))
    continue; %we can't simply average voxels for ipsi and contra conditions
  end
  bvf.(fields{i})=(lvf.(fields{i}) + flipThis(rvf.(fields{i})))/2;
end
%% a hack to flip matrix, should be same as flip in Matlab 2014b, but in earlier versions flip is not available, hence this hack
function y=flipThis(x)
dims=size(x);
if length(dims)>2
  error('this function is not implemented for matrix with more than 2 dims');
end
if dims(1)==1
  y=flipdim(x,2);
else
  y=flipdim(x,1);
end
%% plot the cross contrast channel result
function plotChannelResultXC(s,roiList,chance,saveName)
sess=strtokr(pwd,'/');
figure;
set(gcf,'name',[sess,':Channel results: ',saveName]);
channelPref=[0 22.5 45 67.5 90 112.5 135 157.5 180]; %this is a hack, should figure this out from input
ctLeftLow=s.lvf.fmLow.chanResp;
ctRightLow=s.rvf.fmLow.chanResp;
ctBothLow=s.bvf.fmLow.chanResp;

%duplicate 1st point as it's the same as the last point
ctLeftLow=[ctLeftLow, ctLeftLow(:,1)]; ctRightLow=[ctRightLow, ctRightLow(:,1)]; ctBothLow=[ctBothLow, ctBothLow(:,1)];

ctLeftHigh=s.lvf.fmHigh.chanResp;
ctRightHigh=s.rvf.fmHigh.chanResp;
ctBothHigh=s.bvf.fmHigh.chanResp;

%duplicate 1st point as it's the same as the last point
ctLeftHigh=[ctLeftHigh, ctLeftHigh(:,1)]; ctRightHigh=[ctRightHigh, ctRightHigh(:,1)]; ctBothHigh=[ctBothHigh, ctBothHigh(:,1)];

subplot(4,3,1);
plot(channelPref,ctLeftLow','o-','linewidth',2);
legend(roiList);
title('Low contrast: left VF');

subplot(4,3,2);
plot(channelPref,ctRightLow','o-','linewidth',2);
legend(roiList);
title('Low contrast: right VF');

subplot(4,3,3);
plot(channelPref,ctBothLow','o-','linewidth',2);
legend({'ipsi','contra'});
title('Low contrast: both VF');

subplot(4,3,4);
plot(channelPref,ctLeftHigh','o-','linewidth',2);
legend(roiList);
title('High contrast: left VF');

subplot(4,3,5);
plot(channelPref,ctRightHigh','o-','linewidth',2);
legend(roiList);
title('High contrast: right VF');

subplot(4,3,6);
plot(channelPref,ctBothHigh','o-','linewidth',2);
legend({'ipsi','contra'});
title('High contrast: both VF');

accu=[s.lvf.fmLow.chanClass, s.lvf.fmHigh.chanClass]';
subplot(4,3,7);
bar(accu);
set(gca,'xticklabel',{'lowC','highC'});
yaxis([0 chance*4]);
hline(chance);
title('channel classification left VF');
% legend('lowC','highC');%roiList);
legend(roiList);

accu=[s.rvf.fmLow.chanClass, s.rvf.fmHigh.chanClass]';
subplot(4,3,8);
bar(accu);
set(gca,'xticklabel',{'lowC','highC'});
yaxis([0 chance*4]);
hline(chance);
title('channel classification right VF');
legend(roiList);

accu=[s.bvf.fmLow.chanClass, s.bvf.fmHigh.chanClass]';
subplot(4,3,9);
bar(accu);
set(gca,'xticklabel',{'lowC','highC'});
yaxis([0 chance*4]);
hline(chance);
legend({'ipsi','contra'});
title('channel classification both VF');

subplot(4,3,10);
plot(s.lvf.stimVals',[s.lvf.fmLow.meanIns',s.lvf.fmHigh.meanIns'],'-o','linewidth',2);
legend(['low:',roiList{1}],['low:',roiList{2}],['high:',roiList{1}],['high:',roiList{2}]);
title('Avg instance: left VF');
xlabel('Orientation');
ylabel('%Signal change');

subplot(4,3,11);
plot(s.rvf.stimVals',[s.rvf.fmLow.meanIns',s.rvf.fmHigh.meanIns'],'-o','linewidth',2);
legend(['low:',roiList{1}],['low:',roiList{2}],['high:',roiList{1}],['high:',roiList{2}]);
title('Avg instance: right VF');
xlabel('Orientation');
ylabel('%Signal change');

subplot(4,3,12);
plot(s.rvf.stimVals',[s.bvf.fmLow.meanIns',s.bvf.fmHigh.meanIns'],'-o','linewidth',2);
legend('low:ipsi','low:contra','high:ipsi','high:contra');
title('Avg instance: Both VF');
xlabel('Orientation');
ylabel('%Signal change');

%% plot the r2 values from forward modeling
function plotChannelR2(s,fieldName,roiList,saveName)
sess=strtokr(pwd,'/');
figure;
set(gcf,'name',[sess,':Channel R2: ',saveName,'--',fieldName]);
channelPref=[0 22.5 45 67.5 90 112.5 135 157.5]; %this is a hack, should figure this out from input
r2LeftLow=s.lvf.fmLow.avgR2.(fieldName);
r2RightLow=s.rvf.fmLow.avgR2.(fieldName);
r2BothLow=s.bvf.fmLow.avgR2.(fieldName);

r2LeftHigh=s.lvf.fmHigh.avgR2.(fieldName);
r2RightHigh=s.rvf.fmHigh.avgR2.(fieldName);
r2BothHigh=s.bvf.fmHigh.avgR2.(fieldName);

subplot(3,9,1);
plot(channelPref,r2LeftLow.clAvg','o-','linewidth',2); hold on
title('Low: lvf-classAvg');
legend(roiList);
subplot(3,9,2);
plot(channelPref,r2RightLow.clAvg','o-','linewidth',2);
title('Low: rvf-classAvg');
subplot(3,9,3);
plot(channelPref,r2BothLow.clAvg','o-','linewidth',2);
legend({'ipsi','contra'});
title('Low: bvf-classAvg');

subplot(3,9,4);
plot(channelPref,r2LeftLow.clAcc','*-','linewidth',2); hold on
title('Low: lvf-classAcc');
legend(roiList);
subplot(3,9,5);
plot(channelPref,r2RightLow.clAcc','o-','linewidth',2);
title('Low: rvf-classAcc');
subplot(3,9,6);
plot(channelPref,r2BothLow.clAcc','o-','linewidth',2);
legend({'ipsi','contra'});
title('Low: bvf-classAcc');

subplot(3,9,7);
plot(channelPref,r2LeftLow.clOne','s-','linewidth',2); 
title('Low: lvf-classOne');
legend(roiList);
subplot(3,9,8);
plot(channelPref,r2RightLow.clOne','o-','linewidth',2);
title('Low: rvf-classOne');
subplot(3,9,9);
plot(channelPref,r2BothLow.clOne','o-','linewidth',2);
legend({'ipsi','contra'});
title('Low: bvf-classOne');

subplot(3,9,10);
plot(channelPref,r2LeftHigh.clAvg','o-','linewidth',2); hold on
title('High: lvf-classAvg');
legend(roiList);
subplot(3,9,11);
plot(channelPref,r2RightHigh.clAvg','o-','linewidth',2);
title('High: rvf-classAvg');
subplot(3,9,12);
plot(channelPref,r2BothHigh.clAvg','o-','linewidth',2);
legend({'ipsi','contra'});
title('High: bvf-classAvg');

subplot(3,9,13);
plot(channelPref,r2LeftHigh.clAcc','*-','linewidth',2); hold on
title('High: lvf-classAcc');
legend(roiList);
subplot(3,9,14);
plot(channelPref,r2RightHigh.clAcc','o-','linewidth',2);
title('High: rvf-classAcc');
subplot(3,9,15);
plot(channelPref,r2BothHigh.clAcc','o-','linewidth',2);
legend({'ipsi','contra'});
title('High: bvf-classAcc');

subplot(3,9,16);
plot(channelPref,r2LeftHigh.clOne','s-','linewidth',2); 
legend(roiList);
title('High: lvf-classOne');
subplot(3,9,17);
plot(channelPref,r2RightHigh.clOne','o-','linewidth',2);
title('High: rvf-classOne');
subplot(3,9,18);
plot(channelPref,r2BothHigh.clOne','o-','linewidth',2);
legend({'ipsi','contra'});
title('High: bvf-classOne');

subplot(3,9,19);
bar(r2LeftLow.overall); 
title('Low: lvf-overall');
subplot(3,9,20);
bar(r2RightLow.overall); 
set(gca,'xticklabel',roiList);
title('Low: rvf-overall');
subplot(3,9,21);
bar(r2BothLow.overall); 
set(gca,'xticklabel',{'ipsi','contra'});
title('Low: bvf-overall');

subplot(3,9,22);
bar(r2LeftHigh.overall); 
title('High: lvf-overall');
subplot(3,9,23);
bar(r2RightHigh.overall); 
set(gca,'xticklabel',roiList);
title('High: rvf-overall');
subplot(3,9,24);
bar(r2BothHigh.overall); 
set(gca,'xticklabel',{'ipsi','contra'});
title('High: bvf-overall');


%% this is the plotting function for the previous leave-one-out result
function plotChannelResult(s,fieldName,roiList,chance,saveName)
sess=strtokr(pwd,'/');
figure;
set(gcf,'name',[sess,':Channel results: ',saveName,'--',fieldName]);
channelPref=[0 22.5 45 67.5 90 112.5 135 157.5 180]; %this is a hack, should figure this out from input
ctLeftLow=s.lvf.fmLow.chanResp.(fieldName);
ctRightLow=s.rvf.fmLow.chanResp.(fieldName);
ctBothLow=s.bvf.fmLow.chanResp.(fieldName);

%duplicate 1st point as it's the same as the last point
ctLeftLow=[ctLeftLow, ctLeftLow(:,1)]; ctRightLow=[ctRightLow, ctRightLow(:,1)]; ctBothLow=[ctBothLow, ctBothLow(:,1)];

ctLeftHigh=s.lvf.fmHigh.chanResp.(fieldName);
ctRightHigh=s.rvf.fmHigh.chanResp.(fieldName);
ctBothHigh=s.bvf.fmHigh.chanResp.(fieldName);

%duplicate 1st point as it's the same as the last point
ctLeftHigh=[ctLeftHigh, ctLeftHigh(:,1)]; ctRightHigh=[ctRightHigh, ctRightHigh(:,1)]; ctBothHigh=[ctBothHigh, ctBothHigh(:,1)];

subplot(4,3,1);
plot(channelPref,ctLeftLow','o-','linewidth',2);
legend(roiList);
title('Low contrast: channel tuning left VF');

subplot(4,3,2);
plot(channelPref,ctRightLow','o-','linewidth',2);
legend(roiList);
title('Low contrast: channel tuning right VF');

subplot(4,3,3);
plot(channelPref,ctBothLow','o-','linewidth',2);
legend({'ipsi','contra'});
title('Low contrast: channel tuning both VF');

subplot(4,3,4);
plot(channelPref,ctLeftHigh','o-','linewidth',2);
legend(roiList);
title('High contrast: channel tuning left VF');

subplot(4,3,5);
plot(channelPref,ctRightHigh','o-','linewidth',2);
legend(roiList);
title('High contrast: channel tuning right VF');

subplot(4,3,6);
plot(channelPref,ctBothHigh','o-','linewidth',2);
legend({'ipsi','contra'});
title('High contrast: channel tuning both VF');

accu=[s.lvf.fmLow.chanClass.(fieldName);s.lvf.fmHigh.chanClass.(fieldName)]';
subplot(4,3,7);
bar(accu);
set(gca,'xticklabel',roiList);
yaxis([0 chance*3]);
hline(chance);
title('channel classification left VF');

accu=[s.rvf.fmLow.chanClass.(fieldName);s.rvf.fmHigh.chanClass.(fieldName)]';
subplot(4,3,8);
bar(accu);
set(gca,'xticklabel',roiList);
yaxis([0 chance*3]);
hline(chance);
title('channel classification right VF');

accu=[s.bvf.fmLow.chanClass.(fieldName);s.bvf.fmHigh.chanClass.(fieldName)]';
subplot(4,3,9);
bar(accu);
set(gca,'xticklabel',{'ipsi','contra'});
yaxis([0 chance*3]);
hline(chance);
legend({'lowC','highC'});
title('channel classification both VF');

subplot(4,3,10);
plot(s.lvf.fmLow.stimVals',[s.lvf.fmLow.meanIns.(fieldName)',s.lvf.fmHigh.meanIns.(fieldName)'],'-o','linewidth',2);
legend(['low:',roiList{1}],['low:',roiList{2}],['high:',roiList{1}],['high:',roiList{2}]);
title('Avg instance: left VF');
xlabel('Orientation');
ylabel('%Signal change');

subplot(4,3,11);
plot(s.rvf.fmLow.stimVals',[s.rvf.fmLow.meanIns.(fieldName)',s.rvf.fmHigh.meanIns.(fieldName)'],'-o','linewidth',2);
legend(['low:',roiList{1}],['low:',roiList{2}],['high:',roiList{1}],['high:',roiList{2}]);
title('Avg instance: right VF');
xlabel('Orientation');
ylabel('%Signal change');

subplot(4,3,12);
plot(s.rvf.fmHigh.stimVals',[s.bvf.fmLow.meanIns.(fieldName)',s.bvf.fmHigh.meanIns.(fieldName)'],'-o','linewidth',2);
legend('low:ipsi','low:contra','high:ipsi','high:contra');
title('Avg instance: Both VF');
xlabel('Orientation');
ylabel('%Signal change');

contraChanResp=[ctBothLow(2,:); ctBothHigh(2,:)];
chanFit=fitTuningWithVM(contraChanResp, channelPref,1);

function fittedVals=fitTuningWithVM(chanResp,channelPref,dispFig)

if ~all(chanResp(:,1)==chanResp(:,end))
  disp('First and last element of channel response different, wrapping around it by appending the first element to the end.');
  chanResp=[chanResp,chanResp(:,1)]; %wrap around first and last value
end

paramInit=[pi,1,0,1]; %reasonable guesses of initial value
paramLB=[0, -inf, -inf, 0];
paramUB=[2*pi, inf, inf, inf];
% [paramFit1,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*2),chanResp(1,:),paramLB,paramUB);
% [paramFit2,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*2),chanResp(2,:),paramLB,paramUB);
[paramFit1,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*2),chanResp(1,:));
[paramFit2,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*2),chanResp(2,:));

xp= channelPref(1):1:channelPref(end)*2;
xpOri=xp/2;

yp1=fit2VM(paramFit1,d2r(xp));
yp1_base0=yp1-min(yp1);
rgHalf1=find(yp1_base0>max(yp1_base0)/2);
% fwhm1=xpOri(rgHalf1(end))-xpOri(rgHalf1(1));
fwhm1=length(rgHalf1)/2;

yp2=fit2VM(paramFit2,d2r(xp));
yp2_base0=yp2-min(yp2);
rgHalf2=find(yp2_base0>max(yp2_base0)/2);
fwhm2=xpOri(rgHalf2(end))-xpOri(rgHalf2(1));

%convert the mu from 2pi to pi, because it's orientation
paramFit1(1)=paramFit1(1)/2;
paramFit2(1)=paramFit2(1)/2;
fittedVals=[[paramFit1, fwhm1];[paramFit2, fwhm2]];

if dispFig %plot the fitted von Mises
  figure;
  subplot(1,2,1);
  plot(channelPref, chanResp(1,:),'o'); hold on
  plot(xpOri,yp1,'r--');
  yaxis([0 0.5]);
  line([xpOri(rgHalf1(1)) xpOri(rgHalf1(1))], [0 yp1(rgHalf1(1))]);
  line([xpOri(rgHalf1(end)) xpOri(rgHalf1(end))], [0 yp1(rgHalf1(end))]);
  line([xpOri(rgHalf1(1)) xpOri(rgHalf1(end))], [yp1(rgHalf1(1)), yp1(rgHalf1(end))]);
  title(['mu=',num2str(r2d(paramFit1(1))), ' base=',num2str(paramFit1(3)), ' amp=',num2str(paramFit1(4)),' fwhm=',num2str(fwhm1)]);

  subplot(1,2,2);
  plot(channelPref, chanResp(2,:),'o'); hold on
  plot(xpOri,yp2,'r--');
  yaxis([0 0.5]);
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


%% plot classification accuracy
function plotAllAccu(s, roiList, chance, saveName, contrast)
% plotClassAccu(s.lvf.pcorr,roiList,[saveName, ':left VF'],chance);
% plotClassAccu(s.rvf.pcorr,roiList,[saveName, ':right VF'],chance);
% plotClassAccu(s.bvf.pcorr,{'ipsi','contral'},[saveName, ':combined VF'],chance);
sess=strtokr(pwd,'/');
figure;
set(gcf,'name',[sess,':Classifier performance ',saveName, ' :', contrast]);
idx=1;
fieldName=['pcorr',contrast];
idx=plotClassAccuWOpca(s.lvf.(fieldName),roiList,'Left VF',chance,idx);
idx=plotClassAccuWOpca(s.rvf.(fieldName),roiList,'Right VF',chance,idx);
idx=plotClassAccuWOpca(s.bvf.(fieldName),{'ipsi','contral'},'Combined VF',chance,idx);

%% plot classification accuracy without pca preprocessing
function idx=plotClassAccuWOpca(pcorr, cond, whichROI, chance,idx)

subplot(3,2,idx);
pfisher=[pcorr.fisher; pcorr.fisherDemean; pcorr.fisherZ]';
bar(pfisher);
yaxis([0 min([2.5*chance 1])]);
hline(chance);
set(gca,'xticklabel',cond);
title(['Fisher LDA: ', whichROI]);

idx=idx+1;
subplot(3,2,idx);
pMahal=[pcorr.Mahal; pcorr.MahalDemean; pcorr.MahalZ]';
bar(pMahal);
yaxis([0 2.5*chance]);
hline(chance);
set(gca,'xticklabel',cond);
title(['Mahal: ',whichROI]);
if idx==2
  legend('plain','demean','zscore')
end
idx=idx+1;

%% plot classifier result
function plotClassAccu(pcorr, cond, whichROI, chance)

figure;
set(gcf,'name',whichROI);
subplot(2,2,1);
pfisher=[pcorr.fisher, pcorr.fisherDemean, pcorr.fisherZ];
bar(pfisher);
yaxis([0 1]);
hline(chance);
set(gca,'xticklabel',cond);
title('Fisher LDA');

subplot(2,2,2);
pMahal=[pcorr.Mahal,pcorr.MahalDemean,pcorr.MahalZ];
bar(pMahal);
yaxis([0 1]);
hline(chance);
legend('plain','demean','zscore')
set(gca,'xticklabel',cond);
title('Mahalanobis');

subplot(2,2,3);
plot(pcorr.fpca','o-');
yaxis([0 1]);
hline(chance);
legend(cond);
xlabel('PCA component');
title('Fisher LDA with PCA');

subplot(2,2,4);
plot(pcorr.mpca','o-');
yaxis([0 1]);
hline(chance);
legend(cond);
xlabel('PCA component');
title('Mahalanobis with PCA');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   cross validation classification and forward modeling  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = leaveOneRunOut(s,roiList,stimVolString,cond,varargin)
%TSL: assuming relevant data is in Concatenation group, first scan, and stimvol is in 1st task, 2nd phase, 1st segment
getArgs(varargin,{'groupName=Concatenation','scanNum=1','taskNum=1','phaseNum=2','segmentNum=1','rankMethod=r2','canonicalType=allfit1'});
v = newView;
if isempty(v),return,end
v = viewSet(v,'curGroup','Concatenation');
v = viewSet(v, 'curScan', scanNum);
concatInfo = viewGet(v,'concatInfo');
numRuns=concatInfo.n;

[s.stimvol s.stimNames] = getStimvol(v,stimVolString,'taskNum',1,'phaseNum',phaseNum,'segmentNum',segmentNum);
s.roi = loadROITSeries(v,roiList,[],[],'straightXform=1'); %straightXform=1 retrieves all voxels in the ROI without size-matching, which could delete useful voxels for MVPA purpose
s.nROI = length(roiList);
if strcmp(rankMethod,'r2')
  v = loadAnalysis(v,'erAnal/erAnalLVF'); 
  s.r2Right = viewGet(v,'overlayData',scanNum);
  v = loadAnalysis(v,'erAnal/erAnalRVF'); 
  s.r2Left = viewGet(v,'overlayData',scanNum);
  if isempty(s.r2Right) || isempty(s.r2Left)
    error('Cannot find the r2 map needed to sort the voxels');
  end
  s.roi= getSortIndex(v,s.roi,s.r2Left,s.r2Right);
end
% save testdata s  %for debug purpose only
% load testdata; disp('&&&&&&&&&&&&&&&&&&This is for debug only, loading a test data set&&&&&&&&&&&&&&&&&&');
if s.rebuild==2
  if strcmp(s.instanceMethod,'mean')
    s.roi = getInstances(v,s.roi,s.stimvol,'blockLen=7','type=mean','n=inf'); %TSL:the first arg, v, is not really used but kept to make it backward compatible (see getCanonical)
  elseif strcmp(s.instanceMethod,'deconv')
    s.roi = getInstances(v,s.roi,s.stimvol,'type=deconv','n=inf','canonicalType',canonicalType);
  end
  return;
end

disp('*********(leaveOneRunOut) Default: now doing leave one run out cross validation**********');
disppercent(-1/numRuns);
for iFold=1:numRuns
  %   split data into train and test set
  clear sTrain sTest;
  [sTrain sTest]=splitRuns(s,iFold);
  if strcmp(s.instanceMethod,'mean')
    sTrain.roi = getInstances(v,sTrain.roi,sTrain.stimvol,'blockLen=7','type=mean','n=inf'); %TSL:the first arg, v, is not really used but kept to make it backward compatible (see getCanonical)
  elseif strcmp(s.instanceMethod,'deconv')
    sTrain.roi = getInstances(v,sTrain.roi,sTrain.stimvol,'type=deconv','n=inf','canonicalType',canonicalType);
  end
  sTest.roi = getInstances(v,sTest.roi,sTest.stimvol,sTrain.roi); %get instances for testing data, same call with training ROIs passed in

  %for the mean method, the instances are normalized time course values with
  %1 added as baseline, so need to subtract them out.
  if strcmp(s.instanceMethod, 'mean')
    sTrain=getOverallMeanIns(sTrain,true);
    sTest=getOverallMeanIns(sTest,true);
  elseif strcmp(s.instanceMethod, 'deconv')
    sTrain=getOverallMeanIns(sTrain,false);
    sTest=getOverallMeanIns(sTest,false);
  end

  %Note here and below we are separating high and low contrast instances;
  %we need to get them at the same time and then separate. Cannot get them
  %separately (Ok for mean but not Ok for deconv)
  sTrainLow=subsetInstances(sTrain, [1:8]);
  sTestLow=subsetInstances(sTest, [1:8]);
  sTrainHigh=subsetInstances(sTrain, [9:16]);
  sTestHigh=subsetInstances(sTest, [9:16]);

  %get the grand mean of the instances, reality check of contrast response
  meanIns=cellfun(@(x) x.instance.meanIns, sTestLow.roi,'UniformOutput',false);
  lowContrIns(iFold,:,:)=[meanIns{1};meanIns{2}];
  meanIns=cellfun(@(x) x.instance.meanIns, sTestHigh.roi,'UniformOutput',false);
  highContrIns(iFold,:,:)=[meanIns{1};meanIns{2}];

  if strcmp(s.mode,'WI')
    sTrainForLow=sTrainLow;
    sTrainForHigh=sTrainHigh;
  elseif strcmp(s.mode, 'BB')
    sTrainForLow=combineLowHigh(sTrainLow,sTrainHigh); %combine high and low contrast
    sTrainForHigh=sTrainForLow;
  end
  if s.doClassification
    pLow(iFold)=doOneRoundClassification(sTrainLow, sTestLow);
    pHigh(iFold)=doOneRoundClassification(sTrainHigh, sTestHigh);
  end
  if s.doForwardModel
    [allRespLow(iFold) allr2Low(iFold) allClassifyLow(iFold) stimValsLow]=doOneRoundForwardModel(sTrainForLow,sTestLow, cond);
    [allRespHigh(iFold) allr2High(iFold) allClassifyHigh(iFold) stimValsHigh]=doOneRoundForwardModel(sTrainForHigh,sTestHigh, cond);
  end
  disppercent(iFold/numRuns,sprintf('Finished %i/%i folds of cross-validation',iFold,numRuns));
end
deleteView(v);

if s.doClassification
  pcorrLow=avgClassAccuAcrossCV(pLow);
  pcorrHigh=avgClassAccuAcrossCV(pHigh); 
  s.pcorrLow=pcorrLow;
  s.pcorrHigh=pcorrHigh;
end

if s.doForwardModel
  [fmLow.chanResp fmLow.chanClass fmLow.avgR2 fmLow.meanIns]=avgChannelAcrossCV(allRespLow, allr2Low, allClassifyLow, lowContrIns, s.nROI);
  fmLow.stimVals=stimValsLow;
  [fmHigh.chanResp fmHigh.chanClass fmHigh.avgR2 fmHigh.meanIns]=avgChannelAcrossCV(allRespHigh, allr2High, allClassifyHigh, highContrIns, s.nROI);
  fmHigh.stimVals=stimValsHigh;
  s.fmLow=fmLow;
  s.fmHigh=fmHigh;  
end

%% This function does one round of FEM given training and test data
function [avgChanResp, r2, classAccu, stimVals]=doOneRoundForwardModel(sTrain,sTest, cond)
allStimVals=stimValsFromNames(sTrain.stimNames,'condName',cond);
stimVals=[allStimVals.(cond)];
[avgChanResp.raw, r2.raw, classAccu.raw]=testChannels(sTest.roi,stimVals, buildChannels(sTrain.roi, stimVals));
% [avgChanResp.z, r2.z, classAccu.z]=testChannels(sTest.roi,stimVals, buildChannels(sTrain.roi, stimVals, 'zscore=1','dispChannels=0'));

%% average channel result across validation folds
function [avgResp, avgClassify, avgR2, avgMeanIns]=avgChannelAcrossCV(allAvgResp, allr2, allClassify, allMeanIns, nroi)
  fields=fieldnames(allAvgResp);
  for i=1:length(fields)
    for j=1:length(allAvgResp)
      thisResp(j,:,:)=allAvgResp(j).(fields{i});
      thisAccu(j,:,:)=allClassify(j).(fields{i});
      thisr2_clAvg(j,:,:)=vertcat(allr2(j).(fields{i}).clAvg);
      thisr2_clAcc(j,:,:)=vertcat(allr2(j).(fields{i}).clAcc);
      thisr2_clOne(j,:,:)=vertcat(allr2(j).(fields{i}).clOneshot);
      for k=1:nroi
        thisr2_voxAvg{k}(j,:)=allr2(j).(fields{i})(k).voxAvg;
        thisr2_voxAcc{k}(j,:)=allr2(j).(fields{i})(k).voxAcc;
        thisr2_voxOne{k}(j,:)=allr2(j).(fields{i})(k).voxOneshot;
      end
      thisr2_overall(j,:)=vertcat(allr2(j).(fields{i}).overall);
    end
    avgResp.(fields{i})=squeeze(mean(thisResp,1));
    avgClassify.(fields{i})=sum(thisAccu(:,:,1))./sum(thisAccu(:,:,2));
    avgMeanIns.(fields{i})=squeeze(mean(allMeanIns,1));
    avgR2.(fields{i}).clAvg=squeeze(mean(thisr2_clAvg,1));
    avgR2.(fields{i}).clAcc=squeeze(mean(thisr2_clAcc,1));
    avgR2.(fields{i}).clOne=squeeze(mean(thisr2_clOne,1));
    avgR2.(fields{i}).voxAvg=cellfun(@(x) squeeze(mean(x)), thisr2_voxAvg,'uniformoutput',false);
    avgR2.(fields{i}).voxAcc=cellfun(@(x) squeeze(mean(x)), thisr2_voxAcc,'uniformoutput',false);
    avgR2.(fields{i}).voxOne=cellfun(@(x) squeeze(mean(x)), thisr2_voxOne,'uniformoutput',false);
    avgR2.(fields{i}).overall=squeeze(mean(thisr2_overall,1));
 end
 
%% perform classification with different classifiers and preprocessing
function p=doOneRoundClassification(sTrain, sTest)

p.fisher=testClassifier(sTest.roi,buildClassifier(sTrain.roi,'type=fisher'));
p.fisherDemean=testClassifier(sTest.roi,buildClassifier(sTrain.roi,'demean=1','type=fisher'));
p.fisherZ=testClassifier(sTest.roi,buildClassifier(sTrain.roi,'zscore=1','type=fisher'));
    
p.Mahal = testClassifier(sTest.roi,buildClassifier(sTrain.roi,'type=mahalanobis'));
p.MahalDemean = testClassifier(sTest.roi,buildClassifier(sTrain.roi,'demean=1','type=mahalanobis'));
p.MahalZ = testClassifier(sTest.roi,buildClassifier(sTrain.roi,'zscore=1','type=mahalanobis'));

%% average classification accuracies across validation folds
function pcorr=avgClassAccuAcrossCV(p)
fields=fieldnames(p);
for i=1:length(fields)
  thisAccu=cat(1,p.(fields{i}));
  pcorr.(fields{i})=mean(thisAccu, 1);
end

%% this uses low contrast data to construct channel first, and then test separately on
%low and high contrast. The test on low is circular, but the key is what
%happens to high contrast data.
function s = crossContrastValidation(s,roiList,stimVolStr,cond,varargin)
%TSL: assuming relevant data is in Concatenation group, first scan, and stimvol is in 1st task, 2nd phase, 1st segment
getArgs(varargin,{'groupName=Concatenation','scanNum=1','taskNum=1','phaseNum=2','segmentNum=1','canonicalType=allfit1'});
v = newView;
if isempty(v),return,end
v = viewSet(v,'curGroup','Concatenation');
v = viewSet(v, 'curScan', scanNum);
concatInfo = viewGet(v,'concatInfo');

[s.stimvol s.stimNames] = getStimvol(v,stimVolStr,'taskNum',1,'phaseNum',phaseNum,'segmentNum',segmentNum);
s.roi = loadROITSeries(v,roiList,[],[],'straightXform=1'); %straightXform=1 retrieves all voxels in the ROI without size-matching, which could delete useful voxels for MVPA purpose
s.nROI = length(roiList);
s.roi = getInstances(v,s.roi,s.stimvol,'blockLen=7','n=inf',['type=',s.instanceMethod],['canonicalType=',canonicalType]); %TSL:the first arg, v, is not really used but kept to make it backward compatible (see getCanonical)

%Note here and below we are separating high and low contrast instances;
%we need to get them at the same time and then separate. Cannot get them
%separately (Ok for mean but not Ok for deconv)
sL=subsetInstances(s, [1:8]); %low contrast is the first eight
sH=subsetInstances(s, [9:16]); %high contrast the last eight

%for the mean method, the instances are normalized time course values with
%1 added as baseline, so need to subtract them out.
if strcmp(s.instanceMethod, 'mean')
  sL=getOverallMeanIns(sL,true);
  sH=getOverallMeanIns(sH,true);
elseif strcmp(s.instanceMethod, 'deconv')
  sL=getOverallMeanIns(sL,false);
  sH=getOverallMeanIns(sH,false);
end

%get the grand mean of the instances, reality check of contrast response
meanIns=cellfun(@(x) x.instance.meanIns, sL.roi,'UniformOutput',false);
lowContrIns=[meanIns{1};meanIns{2}];
meanIns=cellfun(@(x) x.instance.meanIns, sH.roi,'UniformOutput',false);
highContrIns=[meanIns{1};meanIns{2}];

temp=stimValsFromNames(sL.stimNames,'condName',cond);
allStimVals=[temp.(cond)];

%build on low contrast, test high contrast
theChannelLow=buildChannels(sL.roi, allStimVals, cond);
[allAvgRespHi thisR2 allClassHi allStimHi allPredHi]=testChannels(sH.roi,allStimVals,theChannelLow);

%build on high contrast, test low contrast
theChannelHigh=buildChannels(sH.roi, allStimVals, cond);
[allAvgRespLow thisR2 allClassLow allStimLow allPredLow]=testChannels(sL.roi,allStimVals,theChannelHigh);

deleteView(v);

  fm.chanResp=allAvgRespLow;
  fm.chanClass=allClassLow(:,1)./allClassLow(:,2);
  fm.actualStim=allStimLow;
  fm.predStim=allPredLow;
  fm.meanIns=lowContrIns;
  s.fmLow=fm;
  
  fm.chanResp=allAvgRespHi;
  fm.chanClass=allClassHi(:,1)./allClassHi(:,2);
  fm.actualStim=allStimHi;
  fm.predStim=allPredHi;
  fm.meanIns=highContrIns;
  s.fmHigh=fm;
  
  s.stimVals=allStimVals;


%% get overall mean and subtract one if asked
function sOut=getOverallMeanIns(sIn, subtractOne)
sOut=sIn;
if subtractOne
  for j=1:length(sIn.roi)
    sOut.roi{j}.instance.instances=cellfun(@(x) minus(x,1),sIn.roi{j}.instance.instances,'UniformOutput',false);
  end
end

for j=1:length(sOut.roi)
  temp=cellfun(@(x) mean(x(:)),sOut.roi{j}.instance.instances,'UniformOutput',false);
  sOut.roi{j}.instance.meanIns=cell2mat(temp);
end  

%% after extracting instances for all trial types, this function put a subset
%of instances (in this case, by contrast) into the s structure. Note esp
%for deconv instance method, extract all and then subset is the correct way
%to do it. For mean instance method, you can just extract a subset.
function  sOut=subsetInstances(sIn, subset)
sOut=sIn;
sOut.stimvol=sIn.stimvol(subset);
sOut.stimNames=sIn.stimNames(subset);

for j=1:length(sIn.roi)
  sOut.roi{j}.instance.instances=sIn.roi{j}.instance.instances(subset);
  sOut.roi{j}.instance.instanceVol=sIn.roi{j}.instance.instanceVol(subset);
  if isfield(sOut.roi{j}.instance,'meanIns')
    sOut.roi{j}.instance.meanIns=sIn.roi{j}.instance.meanIns(subset);
  end
end  

%% combine two instance sets, basicaly used for combining low and high contrast
function sC=combineLowHigh(sL, sH)
sC=sL;
sC=rmfield(sC,'stimvol');
for i=1:length(sL.stimvol)
  sC.stimvol{i}=[sL.stimvol{i}, sH.stimvol{i}];
end
for j=1:length(sC.roi)
  sC.roi{j}.instance.instances=cellfun(@(x,y) cat(1,x,y), sL.roi{j}.instance.instances, sH.roi{j}.instance.instances,'UniformOutput',false);
  sC.roi{j}.instance.instanceVol=cellfun(@(x,y) cat(2,x,y), sL.roi{j}.instance.instanceVol, sH.roi{j}.instance.instanceVol,'UniformOutput',false);
end  

%%