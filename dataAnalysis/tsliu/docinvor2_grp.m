function  docinvor2_grp(tt, conf)
%group analysis for cinvor2
% docinvor2_grp(0, 3); %within-contrast, decon1g instance: current choice for the manuscript May 25, 2016

subjList={'s00520140704/', 's00720150319/', 's01720140718/', 's01920150212/', 's02120150325/', 's02920150330/'};
nsubj=length(subjList);

switch conf %choice of instance method
  case 1
    analName='meanIns'; %'meanIns-1'; %
  case 2
    analName='deconIns';
  case 3
    analName='decon1gIns';
  otherwise
    error('second input must be 1,2,3');
end

switch tt %chocie of cross-validation
  case 0
    disp(['%%%%%%%%%Averaging within-contrast channel tuning using ', analName,'%%%%%%%%%'])
    avgSubjChannel(subjList, nsubj, analName);
  case 1 
    disp(['%%%%%%%%%Across contrast high->low and low->high, channel tuning with ', analName,'%%%%%%%%%']) 
    analName=[analName,'XC'];
    avgXCanal(subjList, nsubj, analName);
  case 2 
    disp(['%%%%%%%%%Across contrast (high+low)->high/low, channel tuning with ', analName,'%%%%%%%%%']) 
    analName=[analName,'BB'];
    avgSubjChannel(subjList, nsubj, analName);  
end

%average forward model results for within contrast, note here we did
%different preprocessing of instances (raw vs. z-score).
function avgSubjChannel(subjList, nsubj, analName)
load([subjList{1},'Anal/',analName]);
roiList={s.lvf.roi{1}.name,s.lvf.roi{2}.name};
fields=fieldnames(s.lvf.fmLow.chanResp);
chance=1/8;
sa=s;  %make a structure sa to store the averaged data
sa=rmfield(sa,{'lvf','rvf','bvf'});

for i=1:nsubj
  clear s;
  thisName=[subjList{i},'Anal/',analName];
  disp(['Loading previously saved analysis: ',thisName]);
  load(thisName);
  thisList={s.lvf.roi{1}.name,s.lvf.roi{2}.name};
  if ~all(cellfun(@(x,y) isequal(x,y), thisList,roiList))
    error('Subjects have different roi names');
  end
  thisFields=fieldnames(s.lvf.fmLow.chanResp);
%   if ~all(cellfun(@(x,y) isequal(x,y), fields,thisFields))
%     error('Subjects have different field names, for different preprocessing');
%   end
  ss(i)=s; %just extract the subject level data
end

for k=1:1 %length(fields) %first field is raw, only looking at raw now
  tf=fields{k};
  sa=averageChannelOnce(ss,sa,tf);
  plotChannelResult(sa,tf,roiList,chance, analName);
  plotChannelR2(sa,tf,roiList, analName);
end


%average forward model resutls for cross contrast, note here we did
%different preprocessing of instances (raw vs. z-score).
function avgXCanal(subjList, nsubj, analName)
load([subjList{1},'Anal/',analName]);
roiList={s.lvf.roi{1}.name,s.lvf.roi{2}.name};
chance=1/8;
sa=s;
sa=rmfield(sa,{'lvf','rvf','bvf'});
sa.lvf.stimVals=s.lvf.stimVals;
sa.rvf.stimVals=s.rvf.stimVals;

for i=1:nsubj
  clear s;
  thisName=[subjList{i},'Anal/',analName];
  disp(['Loading previously saved analysis: ',thisName]);
  load(thisName);
  thisList={s.lvf.roi{1}.name,s.lvf.roi{2}.name};
  if ~all(cellfun(@(x,y) isequal(x,y), thisList,roiList))
    error('Subjects have different roi names');
  end
  
  allChanResp_lvfLow(i,:,:)=s.lvf.fmLow.chanResp;
  allChanResp_lvfHigh(i,:,:)=s.lvf.fmHigh.chanResp;
  allChanResp_rvfLow(i,:,:)=s.rvf.fmLow.chanResp;
  allChanResp_rvfHigh(i,:,:)=s.rvf.fmHigh.chanResp;
  allChanResp_bvfLow(i,:,:)=s.bvf.fmLow.chanResp;
  allChanResp_bvfHigh(i,:,:)=s.bvf.fmHigh.chanResp;

  allChanClass_lvfLow(i,:)=s.lvf.fmLow.chanClass';
  allChanClass_lvfHigh(i,:)=s.lvf.fmHigh.chanClass';
  allChanClass_rvfLow(i,:)=s.rvf.fmLow.chanClass';
  allChanClass_rvfHigh(i,:)=s.rvf.fmHigh.chanClass';
  allChanClass_bvfLow(i,:)=s.bvf.fmLow.chanClass';
  allChanClass_bvfHigh(i,:)=s.bvf.fmHigh.chanClass';

  allMeanIns_lvfLow(i,:,:)=s.lvf.fmLow.meanIns;
  allMeanIns_lvfHigh(i,:,:)=s.lvf.fmHigh.meanIns;
  allMeanIns_rvfLow(i,:,:)=s.rvf.fmLow.meanIns;
  allMeanIns_rvfHigh(i,:,:)=s.rvf.fmHigh.meanIns;
  allMeanIns_bvfLow(i,:,:)=s.bvf.fmLow.meanIns;
  allMeanIns_bvfHigh(i,:,:)=s.bvf.fmHigh.meanIns;
end
[sa.lvf.fmLow.chanResp, sa.lvf.fmLow.chanRespErr]=getMeanErr(allChanResp_lvfLow);
[sa.rvf.fmLow.chanResp, sa.rvf.fmLow.chanRespErr]=getMeanErr(allChanResp_rvfLow);
[sa.bvf.fmLow.chanResp, sa.bvf.fmLow.chanRespErr]=getMeanErr(allChanResp_bvfLow);

[sa.lvf.fmHigh.chanResp, sa.lvf.fmHigh.chanRespErr]=getMeanErr(allChanResp_lvfHigh);
[sa.rvf.fmHigh.chanResp, sa.rvf.fmHigh.chanRespErr]=getMeanErr(allChanResp_rvfHigh);
[sa.bvf.fmHigh.chanResp, sa.bvf.fmHigh.chanRespErr]=getMeanErr(allChanResp_bvfHigh);
 
[sa.lvf.fmLow.chanClass, sa.lvf.fmLow.chanClassErr]=getMeanErr(allChanClass_lvfLow);
[sa.rvf.fmLow.chanClass, sa.rvf.fmLow.chanClassErr]=getMeanErr(allChanClass_rvfLow);
[sa.bvf.fmLow.chanClass, sa.bvf.fmLow.chanClassErr]=getMeanErr(allChanClass_bvfLow);

[sa.lvf.fmHigh.chanClass, sa.lvf.fmHigh.chanClassErr]=getMeanErr(allChanClass_lvfHigh);
[sa.rvf.fmHigh.chanClass, sa.rvf.fmHigh.chanClassErr]=getMeanErr(allChanClass_rvfHigh);
[sa.bvf.fmHigh.chanClass, sa.bvf.fmHigh.chanClassErr]=getMeanErr(allChanClass_bvfHigh);
 
[sa.lvf.fmLow.meanIns, sa.lvf.fmLow.meanInsErr]=getMeanErr(allMeanIns_lvfLow);
[sa.rvf.fmLow.meanIns, sa.rvf.fmLow.meanInsErr]=getMeanErr(allMeanIns_rvfLow);
[sa.bvf.fmLow.meanIns, sa.bvf.fmLow.meanInsErr]=getMeanErr(allMeanIns_bvfLow);

[sa.lvf.fmHigh.meanIns, sa.lvf.fmHigh.meanInsErr]=getMeanErr(allMeanIns_lvfHigh);
[sa.rvf.fmHigh.meanIns, sa.rvf.fmHigh.meanInsErr]=getMeanErr(allMeanIns_rvfHigh);
[sa.bvf.fmHigh.meanIns, sa.bvf.fmHigh.meanInsErr]=getMeanErr(allMeanIns_bvfHigh);
 
plotChannelResultXC(sa,roiList,chance,analName);


%%average channel response across subjects for within-contrast FM
function sa=averageChannelOnce(ss,sa,tf)
nsubj=length(ss);
sa.stimVals=ss(1).lvf.fmLow.stimVals;
sa.channelPref=[0 22.5 45 67.5 90 112.5 135 157.5 180]; %This is a hack, should maybe store this value from buildChannel --TSL
for i=1:nsubj
  if ~all(sa.stimVals==ss(i).lvf.fmLow.stimVals | sa.stimVals==ss(i).rvf.fmLow.stimVals | ...
      sa.stimVals==ss(i).lvf.fmHigh.stimVals | sa.stimVals==ss(i).rvf.fmHigh.stimVals)
    error('stimvals not the same across subjects and conditions, something is wrong');
  end
  allChanResp_lvfLow(i,:,:)=ss(i).lvf.fmLow.chanResp.(tf);
  allChanResp_lvfHigh(i,:,:)=ss(i).lvf.fmHigh.chanResp.(tf);
  allChanResp_rvfLow(i,:,:)=ss(i).rvf.fmLow.chanResp.(tf);
  allChanResp_rvfHigh(i,:,:)=ss(i).rvf.fmHigh.chanResp.(tf);
  allChanResp_bvfLow(i,:,:)=ss(i).bvf.fmLow.chanResp.(tf);
  allChanResp_bvfHigh(i,:,:)=ss(i).bvf.fmHigh.chanResp.(tf);
  allChanFit_bvfLow(i,:,:)=fitTuningWithVM(ss(i).bvf.fmLow.chanResp.(tf), sa.channelPref,0);
  allChanFit_bvfHigh(i,:,:)=fitTuningWithVM(ss(i).bvf.fmHigh.chanResp.(tf), sa.channelPref,0);
  
  allChanClass_lvfLow(i,:)=ss(i).lvf.fmLow.chanClass.(tf);
  allChanClass_lvfHigh(i,:)=ss(i).lvf.fmHigh.chanClass.(tf);
  allChanClass_rvfLow(i,:)=ss(i).rvf.fmLow.chanClass.(tf);
  allChanClass_rvfHigh(i,:)=ss(i).rvf.fmHigh.chanClass.(tf);
  allChanClass_bvfLow(i,:)=ss(i).bvf.fmLow.chanClass.(tf);
  allChanClass_bvfHigh(i,:)=ss(i).bvf.fmHigh.chanClass.(tf);
  
  allR2_lvfLow(i)=ss(i).lvf.fmLow.avgR2.(tf);
  allR2_lvfHigh(i)=ss(i).lvf.fmHigh.avgR2.(tf);
  allR2_rvfLow(i)=ss(i).rvf.fmLow.avgR2.(tf);
  allR2_rvfHigh(i)=ss(i).rvf.fmHigh.avgR2.(tf);
  allR2_bvfLow(i)=ss(i).bvf.fmLow.avgR2.(tf);
  allR2_bvfHigh(i)=ss(i).bvf.fmHigh.avgR2.(tf);

  allMeanIns_lvfLow(i,:,:)=ss(i).lvf.fmLow.meanIns.(tf);
  allMeanIns_lvfHigh(i,:,:)=ss(i).lvf.fmHigh.meanIns.(tf);
  allMeanIns_rvfLow(i,:,:)=ss(i).rvf.fmLow.meanIns.(tf);
  allMeanIns_rvfHigh(i,:,:)=ss(i).rvf.fmHigh.meanIns.(tf);
  allMeanIns_bvfLow(i,:,:)=ss(i).bvf.fmLow.meanIns.(tf);
  allMeanIns_bvfHigh(i,:,:)=ss(i).bvf.fmHigh.meanIns.(tf);
end
[sa.lvf.fmLow.chanResp.(tf), sa.lvf.fmLow.chanRespErr.(tf)]=getMeanErr(allChanResp_lvfLow);
[sa.rvf.fmLow.chanResp.(tf), sa.rvf.fmLow.chanRespErr.(tf)]=getMeanErr(allChanResp_rvfLow);
[sa.bvf.fmLow.chanResp.(tf), sa.bvf.fmLow.chanRespErr.(tf)]=getMeanErr(allChanResp_bvfLow);
[sa.lvf.fmHigh.chanResp.(tf), sa.lvf.fmHigh.chanRespErr.(tf)]=getMeanErr(allChanResp_lvfHigh);
[sa.rvf.fmHigh.chanResp.(tf), sa.rvf.fmHigh.chanRespErr.(tf)]=getMeanErr(allChanResp_rvfHigh);
[sa.bvf.fmHigh.chanResp.(tf), sa.bvf.fmHigh.chanRespErr.(tf)]=getMeanErr(allChanResp_bvfHigh);
 
[sa.lvf.fmLow.chanClass.(tf), sa.lvf.fmLow.chanClassErr.(tf)]=getMeanErr(allChanClass_lvfLow);
[sa.rvf.fmLow.chanClass.(tf), sa.rvf.fmLow.chanClassErr.(tf)]=getMeanErr(allChanClass_rvfLow);
[sa.bvf.fmLow.chanClass.(tf), sa.bvf.fmLow.chanClassErr.(tf)]=getMeanErr(allChanClass_bvfLow);
[sa.lvf.fmHigh.chanClass.(tf), sa.lvf.fmHigh.chanClassErr.(tf)]=getMeanErr(allChanClass_lvfHigh);
[sa.rvf.fmHigh.chanClass.(tf), sa.rvf.fmHigh.chanClassErr.(tf)]=getMeanErr(allChanClass_rvfHigh);
[sa.bvf.fmHigh.chanClass.(tf), sa.bvf.fmHigh.chanClassErr.(tf)]=getMeanErr(allChanClass_bvfHigh);

[sa.lvf.fmLow.meanIns.(tf), sa.lvf.fmLow.meanInsErr.(tf)]=getMeanErr(allMeanIns_lvfLow);
[sa.rvf.fmLow.meanIns.(tf), sa.rvf.fmLow.meanInsErr.(tf)]=getMeanErr(allMeanIns_rvfLow);
[sa.bvf.fmLow.meanIns.(tf), sa.bvf.fmLow.meanInsErr.(tf)]=getMeanErr(allMeanIns_bvfLow);
[sa.lvf.fmHigh.meanIns.(tf), sa.lvf.fmHigh.meanInsErr.(tf)]=getMeanErr(allMeanIns_lvfHigh);
[sa.rvf.fmHigh.meanIns.(tf), sa.rvf.fmHigh.meanInsErr.(tf)]=getMeanErr(allMeanIns_rvfHigh);
[sa.bvf.fmHigh.meanIns.(tf), sa.bvf.fmHigh.meanInsErr.(tf)]=getMeanErr(allMeanIns_bvfHigh);

[sa.bvf.fmLow.chanFit.(tf), sa.bvf.fmLow.chanFitErr.(tf)]=getMeanWiErr(allChanFit_bvfLow,2);
[sa.bvf.fmHigh.chanFit.(tf), sa.bvf.fmHigh.chanFitErr.(tf)]=getMeanWiErr(allChanFit_bvfHigh,2);

contraIdx=2;
mus=[allChanFit_bvfLow(:,contraIdx,1),allChanFit_bvfHigh(:,contraIdx,1)];
bases=[allChanFit_bvfLow(:,contraIdx,3),allChanFit_bvfHigh(:,contraIdx,3)];
amps=[allChanFit_bvfLow(:,contraIdx,4),allChanFit_bvfHigh(:,contraIdx,4)];
fwhms=[allChanFit_bvfLow(:,contraIdx,5),allChanFit_bvfHigh(:,contraIdx,5)];

[h_mu,p_mu,ci_mu,stats_mu]=ttest(mus(:,1), mus(:,2));
[h_base,p_base,ci_base,stats_base]=ttest(bases(:,1), bases(:,2));
[h_amp,p_amp,ci_amp,stats_amp]=ttest(amps(:,1), amps(:,2));
[h_fwhm,p_fwhm,ci_fwhm,stats_fwhm]=ttest(fwhms(:,1), fwhms(:,2));
disp('Testing channel fit resutls between low vs. high contrast');
disp(['p_val of mus, baseline, amplitude, fwhm=', num2str([p_mu, p_base, p_amp, p_fwhm])]);

sa.lvf.fmLow.avgR2.(tf)=getMeanR2(allR2_lvfLow, nsubj);
sa.rvf.fmLow.avgR2.(tf)=getMeanR2(allR2_rvfLow, nsubj);
sa.bvf.fmLow.avgR2.(tf)=getMeanR2(allR2_bvfLow, nsubj);
sa.lvf.fmHigh.avgR2.(tf)=getMeanR2(allR2_lvfHigh, nsubj);
sa.rvf.fmHigh.avgR2.(tf)=getMeanR2(allR2_rvfHigh, nsubj);
sa.bvf.fmHigh.avgR2.(tf)=getMeanR2(allR2_bvfHigh, nsubj);

allr2_bvfLow_overall=cat(1,allR2_bvfLow.overall);
allr2_bvfHigh_overall=cat(1,allR2_bvfHigh.overall);
[h,p,ci,stats]=ttest(allr2_bvfLow_overall(:,2), allr2_bvfHigh_overall(:,2));
disp(['Testing diff in contralateral r2 across subjects, using overall values. p=',num2str(p)]);

function fittedVals=fitTuningWithVM(chanResp,channelPref,dispFig)

if ~all(chanResp(:,1)==chanResp(:,end))
  disp('First and last element of channel response different, wrapping around it by appending the first element to the end.');
  chanResp=[chanResp,chanResp(:,1)]; %wrap around first and last value
end

paramInit=[pi,1,0,1]; %reasonable guesses of initial value
% paramLB=[0, -inf, -inf, 0];
% paramUB=[2*pi, inf, inf, inf];
% [paramFit1,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*2),chanResp(1,:),paramLB,paramUB);
% [paramFit2,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*2),chanResp(2,:),paramLB,paramUB);
[paramFit1,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*2),chanResp(1,:));
[paramFit2,resnorm,residual,exitflag]=lsqcurvefit(@fit2VM,paramInit,d2r(channelPref*2),chanResp(2,:));

xp= channelPref(1):1:channelPref(end)*2;
xpOri=xp/2;

yp1=fit2VM(paramFit1,d2r(xp));
yp1_base0=yp1-min(yp1);
rgHalf1=find(yp1_base0>max(yp1_base0)/2);
if paramFit1(4)>0
  fwhm1=xpOri(rgHalf1(end))-xpOri(rgHalf1(1));
else
  fwhm1=length(rgHalf1)/2;
end

yp2=fit2VM(paramFit2,d2r(xp));
yp2_base0=yp2-min(yp2);
rgHalf2=find(yp2_base0>max(yp2_base0)/2);
if paramFit2(4)>0
  fwhm2=xpOri(rgHalf2(end))-xpOri(rgHalf2(1));
else
  fwhm2=length(rgHalf2)/2;
end
%convert the mu from 2pi to pi, because it's orientation
paramFit1(1)=paramFit1(1)/2;
paramFit2(1)=paramFit2(1)/2;
fittedVals=[[paramFit1, fwhm1];[paramFit2, fwhm2]];

if dispFig %plot the fitted von Mises
  figure;
  subplot(1,3,1);
  plot(channelPref, chanResp(1,:),'o-'); hold on
  plot(xpOri,yp1,'r--');
  yaxis([0 0.5]);
  line([xpOri(rgHalf1(1)) xpOri(rgHalf1(1))], [0 yp1(rgHalf1(1))]);
  line([xpOri(rgHalf1(end)) xpOri(rgHalf1(end))], [0 yp1(rgHalf1(end))]);
  line([xpOri(rgHalf1(1)) xpOri(rgHalf1(end))], [yp1(rgHalf1(1)), yp1(rgHalf1(end))]);
  title(['mu=',num2str(r2d(paramFit1(1))), ' base=',num2str(paramFit1(3)), ' amp=',num2str(paramFit1(4)),' fwhm=',num2str(fwhm1)]);

  subplot(1,3,2);
  plot(channelPref, chanResp(2,:),'o-'); hold on
  plot(xpOri,yp2,'r--');
  yaxis([0 0.5]);
  line([xpOri(rgHalf2(1)) xpOri(rgHalf2(1))], [0 yp2(rgHalf2(1))]);
  line([xpOri(rgHalf2(end)) xpOri(rgHalf2(end))], [0 yp2(rgHalf2(end))]);
  line([xpOri(rgHalf2(1)) xpOri(rgHalf2(end))], [yp2(rgHalf2(1)), yp2(rgHalf2(end))]);
  title(['mu=',num2str(r2d(paramFit2(1))), ' base=',num2str(paramFit2(3)), ' amp=',num2str(paramFit2(4)),' fwhm=',num2str(fwhm2)]);
  
  hc=[0.8 0.3 0];
  lc=[0 0.4 0.8];
  subplot(1,3,3);
  p1=plot(channelPref, chanResp(1,:),'o'); hold on
  p2=plot(xpOri,yp1,'-','linewidth',2);
  set(p1,'Color',lc);
  set(p2,'Color',lc);
  p3=plot(channelPref, chanResp(2,:),'o'); hold on
  p4=plot(xpOri,yp2,'-','linewidth',2);
  set(p3,'Color',hc);
  set(p4,'Color',hc);
  axis([-10 190 -0.05 0.45]);
  legend({'low','high'});
%   title('Low contrast');
  xlabel('Orientation (deg)');
  ylabel('Channel response (arb)');
  box off;

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

function avgR2=getMeanR2(allR2,nsubj)
fields=fieldnames(allR2);
for i=1:length(fields)
  thisR2=[];
  if strfind(fields{i},'vox')
    continue; %note we're not averaging voxel r2, not clear how to do it for individual voxels    
  else
    for j=1:nsubj
      thisR2(j,:,:)=allR2(j).(fields{i});
    end
    [avgR2.(fields{i}), avgR2.([fields{i},'Err'])]=getMeanErr(thisR2);
  end
end


function plotChannelResultXC(s,roiList,chance,saveName)
figure;
set(gcf,'name',['Average XC results: ',saveName]);
channelPref=[0 22.5 45 67.5 90 112.5 135 157.5 180]; %this is a hack, should figure this out from input
ctLeftLow=s.lvf.fmLow.chanResp; ctLeftLowErr=s.lvf.fmLow.chanRespErr;
ctRightLow=s.rvf.fmLow.chanResp; ctRightLowErr=s.rvf.fmLow.chanRespErr;
ctBothLow=s.bvf.fmLow.chanResp; ctBothLowErr=s.bvf.fmLow.chanRespErr;

%duplicate 1st point as it's the same as the last point
ctLeftLow=[ctLeftLow, ctLeftLow(:,1)]; ctRightLow=[ctRightLow, ctRightLow(:,1)]; ctBothLow=[ctBothLow, ctBothLow(:,1)];
ctLeftLowErr=[ctLeftLowErr, ctLeftLowErr(:,1)]; ctRightLowErr=[ctRightLowErr, ctRightLowErr(:,1)]; ctBothLowErr=[ctBothLowErr, ctBothLowErr(:,1)];

ctLeftHigh=s.lvf.fmHigh.chanResp; ctLeftHighErr=s.lvf.fmHigh.chanRespErr;
ctRightHigh=s.rvf.fmHigh.chanResp; ctRightHighErr=s.rvf.fmHigh.chanRespErr;
ctBothHigh=s.bvf.fmHigh.chanResp; ctBothHighErr=s.bvf.fmHigh.chanRespErr;

%duplicate 1st point as it's the same as the last point
ctLeftHigh=[ctLeftHigh, ctLeftHigh(:,1)]; ctRightHigh=[ctRightHigh, ctRightHigh(:,1)]; ctBothHigh=[ctBothHigh, ctBothHigh(:,1)];
ctLeftHighErr=[ctLeftHighErr, ctLeftHighErr(:,1)]; ctRightHighErr=[ctRightHighErr, ctRightHighErr(:,1)]; ctBothHighErr=[ctBothHighErr, ctBothHighErr(:,1)];

xval=repmat(channelPref',1,2);

subplot(4,3,1);
errorbar(xval,ctLeftLow',ctLeftLowErr','o-','linewidth',2);
xaxis([-10 190]);
legend(roiList);
title('Low contrast: left VF');

subplot(4,3,2);
errorbar(xval,ctRightLow',ctRightLowErr','o-','linewidth',2);
xaxis([-10 190]);
legend(roiList);
title('Low contrast: right VF');

subplot(4,3,3);
errorbar(xval,ctBothLow',ctBothLowErr', 'o-','linewidth',2);
xaxis([-10 190]);
legend({'ipsi','contra'});
title('Low contrast: both VF');

subplot(4,3,4);
errorbar(xval,ctLeftHigh',ctLeftHighErr','o-','linewidth',2);
xaxis([-10 190]);
legend(roiList);
title('High contrast: left VF');

subplot(4,3,5);
errorbar(xval,ctRightHigh',ctRightHighErr', 'o-','linewidth',2);
xaxis([-10 190]);
legend(roiList);
title('High contrast: right VF');

subplot(4,3,6);
errorbar(xval,ctBothHigh',ctBothHighErr', 'o-','linewidth',2);
xaxis([-10 190]);
legend({'ipsi','contra'});
title('High contrast: both VF');

accu=[s.lvf.fmLow.chanClass; s.lvf.fmHigh.chanClass]';
subplot(4,3,7);
bar(accu);
set(gca,'xticklabel',{'lowC','highC'});
yaxis([0 chance*4]);
hline(chance);
title('channel classification left VF');
legend(roiList);

accu=[s.rvf.fmLow.chanClass; s.rvf.fmHigh.chanClass]';
subplot(4,3,8);
bar(accu);
set(gca,'xticklabel',{'lowC','highC'});
yaxis([0 chance*4]);
hline(chance);
title('channel classification right VF');
legend(roiList);

accu=[s.bvf.fmLow.chanClass; s.bvf.fmHigh.chanClass]';
subplot(4,3,9);
bar(accu);
set(gca,'xticklabel',{'lowC','highC'});
yaxis([0 chance*4]);
hline(chance);
legend({'ipsi','contra'});
title('channel classification both VF');

xval=repmat(s.lvf.stimVals',1,4);
subplot(4,3,10); 
errorbar(xval,[s.lvf.fmLow.meanIns',s.lvf.fmHigh.meanIns'],[s.lvf.fmLow.meanInsErr',s.lvf.fmHigh.meanInsErr'],'-o','linewidth',2);
xaxis([-10 190]);
legend(['low:',roiList{1}],['low:',roiList{2}],['high:',roiList{1}],['high:',roiList{2}]);
title('Avg instance: left VF');

subplot(4,3,11);
errorbar(xval,[s.rvf.fmLow.meanIns',s.rvf.fmHigh.meanIns'],[s.rvf.fmLow.meanInsErr',s.rvf.fmHigh.meanInsErr'],'-o','linewidth',2);
xaxis([-10 190]);
legend(['low:',roiList{1}],['low:',roiList{2}],['high:',roiList{1}],['high:',roiList{2}]);
title('Avg instance: right VF');

subplot(4,3,12);
errorbar(xval,[s.bvf.fmLow.meanIns',s.bvf.fmHigh.meanIns'],[s.bvf.fmLow.meanInsErr',s.bvf.fmHigh.meanInsErr'],'-o','linewidth',2);
legend('low:ipsi','low:contra','high:ipsi','high:contra');
title('Avg instance: Both VF');
xlabel('Orientation');
ylabel('%Signal change');


%this is the plotting function for the previous leave-one-out result
function plotChannelResult(s,fieldName,roiList,chance,saveName)
sess=strtokr(pwd,'/');
figure;
set(gcf,'name',[sess,': ',saveName,'--',fieldName]);
channelPref=s.channelPref;
ctLeftLow=s.lvf.fmLow.chanResp.(fieldName); ctLeftLowErr=s.lvf.fmLow.chanRespErr.(fieldName);
ctRightLow=s.rvf.fmLow.chanResp.(fieldName); ctRightLowErr=s.rvf.fmLow.chanRespErr.(fieldName);
ctBothLow=s.bvf.fmLow.chanResp.(fieldName); ctBothLowErr=s.bvf.fmLow.chanRespErr.(fieldName);

%duplicate 1st point as it's the same as the last point
ctLeftLow=[ctLeftLow, ctLeftLow(:,1)]; ctRightLow=[ctRightLow, ctRightLow(:,1)]; ctBothLow=[ctBothLow, ctBothLow(:,1)];
ctLeftLowErr=[ctLeftLowErr, ctLeftLowErr(:,1)]; ctRightLowErr=[ctRightLowErr, ctRightLowErr(:,1)]; ctBothLowErr=[ctBothLowErr, ctBothLowErr(:,1)];

ctLeftHigh=s.lvf.fmHigh.chanResp.(fieldName); ctLeftHighErr=s.lvf.fmHigh.chanRespErr.(fieldName);
ctRightHigh=s.rvf.fmHigh.chanResp.(fieldName); ctRightHighErr=s.rvf.fmHigh.chanRespErr.(fieldName);
ctBothHigh=s.bvf.fmHigh.chanResp.(fieldName); ctBothHighErr=s.bvf.fmHigh.chanRespErr.(fieldName);

%duplicate 1st point as it's the same as the last point
ctLeftHigh=[ctLeftHigh, ctLeftHigh(:,1)]; ctRightHigh=[ctRightHigh, ctRightHigh(:,1)]; ctBothHigh=[ctBothHigh, ctBothHigh(:,1)];
ctLeftHighErr=[ctLeftHighErr, ctLeftHighErr(:,1)]; ctRightHighErr=[ctRightHighErr, ctRightHighErr(:,1)]; ctBothHighErr=[ctBothHighErr, ctBothHighErr(:,1)];

xval=repmat(channelPref',1,2);
subplot(4,3,1);
errorbar(xval,ctLeftLow',ctLeftLowErr','o-','linewidth',2);
xaxis([-10 190]);
legend(roiList);
title('Low contrast: left VF');

subplot(4,3,2);
errorbar(xval,ctRightLow',ctRightLowErr', 'o-','linewidth',2);
xaxis([-10 190]);
legend(roiList);
title('Low contrast: right VF');

subplot(4,3,3);
errorbar(xval,ctBothLow',ctBothLowErr', 'o-','linewidth',2); hold on
axis([-10 190 0 0.5]);
legend({'ipsi','contra'});
title('Low contrast: both VF');

subplot(4,3,4);
errorbar(xval,ctLeftHigh',ctLeftHighErr', 'o-','linewidth',2);
xaxis([-10 190]);
legend(roiList);
title('High contrast: left VF');

subplot(4,3,5);
errorbar(xval,ctRightHigh',ctRightHighErr', 'o-','linewidth',2);
xaxis([-10 190]);
legend(roiList);
title('High contrast: right VF');

subplot(4,3,6);
errorbar(xval,ctBothHigh',ctBothHighErr', 'o-','linewidth',2); hold on
axis([-10 190 0 0.5]);
legend({'ipsi','contra'});
title('High contrast: both VF');

accu=[s.lvf.fmLow.chanClass.(fieldName); s.lvf.fmHigh.chanClass.(fieldName)]';
accuErr=[s.lvf.fmLow.chanClassErr.(fieldName); s.lvf.fmHigh.chanClassErr.(fieldName)]';

subplot(4,3,7);
bars(accu,accuErr,roiList);
% mybar(accu,'yerror',accuErr,'groupLabels',roiList,'withinGroupLabels',{'low','high'});
yaxis([0 chance*3]);
hline(chance);
title('channel classification left VF');

accu=[s.rvf.fmLow.chanClass.(fieldName);s.rvf.fmHigh.chanClass.(fieldName)]';
accuErr=[s.rvf.fmLow.chanClassErr.(fieldName);s.rvf.fmHigh.chanClassErr.(fieldName)]';
subplot(4,3,8);
bars(accu,accuErr,roiList);
yaxis([0 chance*3]);
hline(chance);
title('channel classification right VF');

accu=[s.bvf.fmLow.chanClass.(fieldName);s.bvf.fmHigh.chanClass.(fieldName)]';
accuErr=[s.bvf.fmLow.chanClassErr.(fieldName);s.bvf.fmHigh.chanClassErr.(fieldName)]';
subplot(4,3,9);
bars(accu,accuErr,{'ipsi','contra'}, char('lowC','highC'));
yaxis([0 chance*3]);
hline(chance);
title('channel classification both VF');

xvalOri=repmat(s.stimVals',1,4);
subplot(4,3,10); 
errorbar(xvalOri,[s.lvf.fmLow.meanIns.(fieldName)',s.lvf.fmHigh.meanIns.(fieldName)'],[s.lvf.fmLow.meanInsErr.(fieldName)',s.lvf.fmHigh.meanInsErr.(fieldName)'],'-o','linewidth',2);
xaxis([-10 190]);
legend(['low:',roiList{1}],['low:',roiList{2}],['high:',roiList{1}],['high:',roiList{2}]);
title('Avg instance: left VF');

subplot(4,3,11);
errorbar(xvalOri,[s.rvf.fmLow.meanIns.(fieldName)',s.rvf.fmHigh.meanIns.(fieldName)'],[s.rvf.fmLow.meanInsErr.(fieldName)',s.rvf.fmHigh.meanInsErr.(fieldName)'],'-o','linewidth',2);
xaxis([-10 190]);
legend(['low:',roiList{1}],['low:',roiList{2}],['high:',roiList{1}],['high:',roiList{2}]);
title('Avg instance: right VF');

subplot(4,3,12);
errorbar(xvalOri,[s.bvf.fmLow.meanIns.(fieldName)',s.bvf.fmHigh.meanIns.(fieldName)'],[s.bvf.fmLow.meanInsErr.(fieldName)',s.bvf.fmHigh.meanInsErr.(fieldName)'],'-o','linewidth',2);
xaxis([-10 190]);
legend('low:ipsi','low:contra','high:ipsi','high:contra');
title('Avg instance: Both VF');
xlabel('Orientation');
ylabel('%Signal change');

%generating figure for manuscript
hc=[0.8 0.3 0];
lc=[0 0.4 0.8];
fig2h=figure;
set(fig2h,'Name','Figure 2')
subplot(2,2,1);cla
h1=errorbar(xval(:,1),ctBothLow(1,:)',ctBothLowErr(1,:)', 'o--','linewidth',2); hold on
h2=errorbar(xval(:,2),ctBothLow(2,:)',ctBothLowErr(2,:)', 'o-','linewidth',2);
set(h1,'Color',lc);
set(h2,'Color',lc);
axis([-10 190 -0.05 0.45]);
legend({'ipsi','contra'});
title('Low contrast');
xlabel('Orientation (deg)');
ylabel('Channel response (arb)');
box off;

subplot(2,2,2);
% errorbar(xval,ctBothHigh',ctBothHighErr', 'o-','linewidth',2); hold on
h3=errorbar(xval(:,1),ctBothHigh(1,:)',ctBothHighErr(1,:)', 'o--','linewidth',2); hold on
h4=errorbar(xval(:,2),ctBothHigh(2,:)',ctBothHighErr(2,:)', 'o-','linewidth',2); 
set(h3,'Color',hc);
set(h4,'Color',hc);
axis([-10 190 -0.05 0.45]);
legend({'ipsi','contra'});
title('High contrast');
xlabel('Orientation (deg)');
ylabel('Channel response (arb)');
box off;

ins=[s.bvf.fmLow.meanIns.(fieldName)',s.bvf.fmHigh.meanIns.(fieldName)'];
insErr=[s.bvf.fmLow.meanInsErr.(fieldName)',s.bvf.fmHigh.meanInsErr.(fieldName)'];
subplot(2,2,3); cla
h1=errorbar(xvalOri(:,1),ins(:,1),insErr(:,1),'--o','linewidth',2); hold on
h2=errorbar(xvalOri(:,2),ins(:,2),insErr(:,2),'-o','linewidth',2);
h3=errorbar(xvalOri(:,3),ins(:,3),insErr(:,3),'--o','linewidth',2);
h4=errorbar(xvalOri(:,4),ins(:,4),insErr(:,4),'-o','linewidth',2);
set(h1,'Color',lc);
set(h2,'Color',lc);
set(h3,'Color',hc);
set(h4,'Color',hc);

axis([-10 190 0.5 2.4]);
legend('low:ipsi','low:contra','high:ipsi','high:contra');
title('Mean BOLD response');
xlabel('Orientation');
ylabel('% Signal change');

% fit the average channel response to a von Mises
avgChanFit=fitTuningWithVM([ctBothLow(2,:);ctBothHigh(2,:)], channelPref,1);
disp('Fitted params of channel resp (LowC/HighC), mean, kappa, baseline, amplitude, FWHM');
disp(num2str(avgChanFit));
  
fig3h=figure;
set(gcf,'name',[sess,': ',saveName,'--',fieldName,':channel Fit:Contralateral ROI']);
lowCFit=s.bvf.fmLow.chanFit.(fieldName)(2,:); %note the second row is the contralateral response
highCFit=s.bvf.fmHigh.chanFit.(fieldName)(2,:);
lowCFitErr=s.bvf.fmLow.chanFitErr.(fieldName)(2,:);
highCFitErr=s.bvf.fmHigh.chanFitErr.(fieldName)(2,:);

subplot(2,2,3);
res=r2d([lowCFit(1);highCFit(1)]);
resErr=r2d([lowCFitErr(1);highCFitErr(1)]);
mybar(res,'yError',resErr,'groupLabels',{'Low','High'},'xLabelText','Mean','yLabelText','deg','yAxisMin=0','yAxisMax=120','groupColors',{lc,hc},'dispValues',0);

subplot(2,2,1);
res=[lowCFit(3);highCFit(3)];
resErr=[lowCFitErr(3);highCFitErr(3)];
mybar(res,'yError',resErr,'groupLabels',{'Low','High'},'xLabelText','Baseline','yLabelText','arb','yAxisMin=0','yAxisMax=0.35','groupColors',{lc,hc},'dispValues',0);

subplot(2,2,2);
res=[lowCFit(4);highCFit(4)];
resErr=[lowCFitErr(4);highCFitErr(4)];
mybar(res,'yError',resErr,'groupLabels',{'Low','High'},'xLabelText','Amplitude','yLabelText','arb','yAxisMin=0','yAxisMax=0.35','groupColors',{lc,hc},'dispValues',0);

subplot(2,2,4);
res=[lowCFit(5);highCFit(5)];
resErr=[lowCFitErr(5);highCFitErr(5)];
mybar(res,'yError',resErr,'groupLabels',{'Low','High'},'xLabelText','FWHM','yLabelText','deg','yAxisMin=0','yAxisMax=120','groupColors',{lc,hc},'dispValues',0);

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
errorbar(repmat(channelPref',1,2),r2LeftLow.clAvg',r2LeftLow.clAvgErr','o-','linewidth',2); 
title('Low: lvf-classAvg'); xaxis([-20 200]);
legend(roiList);
subplot(3,9,2);
errorbar(repmat(channelPref',1,2),r2RightLow.clAvg',r2RightLow.clAvgErr','o-','linewidth',2); 
title('Low: rvf-classAvg'); xaxis([-20 200]);
subplot(3,9,3);
errorbar(repmat(channelPref',1,2),r2BothLow.clAvg',r2BothLow.clAvgErr','o-','linewidth',2); 
legend({'ipsi','contra'}); xaxis([-20 200]);
title('Low: bvf-classAvg');

subplot(3,9,4);
errorbar(repmat(channelPref',1,2),r2LeftLow.clAcc',r2LeftLow.clAccErr','o-','linewidth',2);
title('Low: lvf-classAcc'); xaxis([-20 200]);
legend(roiList);
subplot(3,9,5);
errorbar(repmat(channelPref',1,2),r2RightLow.clAcc',r2RightLow.clAccErr','o-','linewidth',2); 
title('Low: rvf-classAcc'); xaxis([-20 200]);
subplot(3,9,6);
errorbar(repmat(channelPref',1,2),r2BothLow.clAcc',r2BothLow.clAccErr','o-','linewidth',2); 
legend({'ipsi','contra'}); xaxis([-20 200]);
title('Low: bvf-classAcc');

subplot(3,9,7);
errorbar(repmat(channelPref',1,2),r2LeftLow.clOne',r2LeftLow.clOneErr','o-','linewidth',2); 
title('Low: lvf-classOne'); xaxis([-20 200]);
legend(roiList);
subplot(3,9,8);
errorbar(repmat(channelPref',1,2),r2RightLow.clOne',r2RightLow.clOneErr','o-','linewidth',2); 
title('Low: rvf-classOne'); xaxis([-20 200]);
subplot(3,9,9);
errorbar(repmat(channelPref',1,2),r2BothLow.clOne',r2BothLow.clOneErr','o-','linewidth',2); 
legend({'ipsi','contra'}); xaxis([-20 200]);
title('Low: bvf-classOne');

subplot(3,9,10);
errorbar(repmat(channelPref',1,2),r2LeftHigh.clAvg',r2LeftHigh.clAvgErr','o-','linewidth',2); 
title('High: lvf-classAvg'); xaxis([-20 200]);
legend(roiList);
subplot(3,9,11);
errorbar(repmat(channelPref',1,2),r2RightHigh.clAvg',r2RightHigh.clAvgErr','o-','linewidth',2); 
title('High: rvf-classAvg'); xaxis([-20 200]);
subplot(3,9,12);
errorbar(repmat(channelPref',1,2),r2BothHigh.clAvg',r2BothHigh.clAvgErr','o-','linewidth',2); 
legend({'ipsi','contra'});
title('High: bvf-classAvg'); xaxis([-20 200]);

subplot(3,9,13);
errorbar(repmat(channelPref',1,2),r2LeftHigh.clAcc',r2LeftHigh.clAccErr','o-','linewidth',2); 
title('High: lvf-classAcc'); xaxis([-20 200]);
legend(roiList);
subplot(3,9,14);
plot(channelPref,r2RightHigh.clAcc','o-','linewidth',2);
errorbar(repmat(channelPref',1,2),r2RightHigh.clAcc',r2RightHigh.clAccErr','o-','linewidth',2); 
title('High: rvf-classAcc'); xaxis([-20 200]);
subplot(3,9,15);
errorbar(repmat(channelPref',1,2),r2BothHigh.clAcc',r2BothHigh.clAccErr','o-','linewidth',2); 
legend({'ipsi','contra'}); xaxis([-20 200]);
title('High: bvf-classAcc');

subplot(3,9,16);
errorbar(repmat(channelPref',1,2),r2LeftHigh.clOne',r2LeftHigh.clOneErr','o-','linewidth',2); 
legend(roiList); xaxis([-20 200]);
title('High: lvf-classOne');
subplot(3,9,17);
errorbar(repmat(channelPref',1,2),r2RightHigh.clOne',r2RightHigh.clOneErr','o-','linewidth',2); 
title('High: rvf-classOne'); xaxis([-20 200]);
subplot(3,9,18);
errorbar(repmat(channelPref',1,2),r2BothHigh.clOne',r2BothHigh.clOneErr','o-','linewidth',2); 
legend({'ipsi','contra'}); xaxis([-20 200]);
title('High: bvf-classOne');

subplot(3,9,19);
mybar(r2LeftLow.overall,'yError',r2LeftLow.overallErr,'groupLabels',roiList);
title('Low: lvf-overall');
subplot(3,9,20);
mybar(r2RightLow.overall,'yError',r2RightLow.overallErr,'groupLabels',roiList);
title('Low: rvf-overall');
subplot(3,9,21);
mybar(r2BothLow.overall,'yError',r2BothLow.overallErr,'groupLabels',{'ipsi','contra'});
title('Low: bvf-overall');

subplot(3,9,22);
mybar(r2LeftHigh.overall,'yError',r2LeftHigh.overallErr,'groupLabels',roiList);
title('High: lvf-overall');
subplot(3,9,23);
mybar(r2RightHigh.overall,'yError',r2RightHigh.overallErr,'groupLabels',roiList);
title('High: rvf-overall');
subplot(3,9,24);
mybar(r2BothHigh.overall,'yError',r2BothHigh.overallErr,'groupLabels',{'ipsi','contra'});
title('High: bvf-overall');

%generating figure for manuscript
hc=[0.8 0.3 0];
lc=[0 0.4 0.8];
fh34=figure;
set(fh34,'Name','Fiure3,panel4');
r2plot=[r2BothLow.overall(2), r2BothHigh.overall(2)];
r2plotErr=[r2BothLow.overallErr(2), r2BothHigh.overallErr(2)];
mybar(r2plot','yError',r2plotErr','groupLabels',{'Low','High'},'yLabelText','r2','groupColors',{lc,hc},'dispValues',0);
