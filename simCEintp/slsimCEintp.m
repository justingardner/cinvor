

%simCEinterp.m
%
% author : steeve laquitaine
%   date : 160831
%purpose : simulate channel encoding reconstruction of 8 stimulus 
%          orientations from 50 response instances of 2 voxels.
%          This tests the case when stimulus orientations do not match the
%          hypothetical channels orientation preferences. Channel responses
%          are interpolated for 360 channels with orientation preferences
%          ranging from 1:1:360 degs.
%          
%
%library required : mgl, mrTools, cinvor, gru, slcircplots
%https://github.com/justingardner/
%https://github.com/steevelaquitaine/


%simulate 50 instances of 2 voxels responses to 36 stimulus orientations
%(the more orientations the better: the more the channel weights are  
%constrained when fitting) when the displayed orientations do not match the 
%orientation preferences of the channels that are used to model the
%voxel responses. Increasing the number of instances doesn't affect the 
%fit because the noise in voxel response has been set to be very low
stimValues = 0:20:350;
numInstances = 100;

%simulate instances for training and train the channel weights
instancesTrain = slsimInst(stimValues,numInstances);
channel = buildChannels(instancesTrain,stimValues);

%simulate test instances and reconstruct orientations from those instances
instancesTest = slsimInst(stimValues,numInstances);
[avgTestResponse r2 classifyCorrTotal stimValVector predStimVal] = testChannels(instancesTest,stimValues,channel,'interpChanResp=1');

%plot linearized reconstructed vs true orientations
%for vizualization and linear fit
sldrawCircStat(predStimVal',stimValVector');

%plot average test channel responses
figure; plot(avgTestResponse)