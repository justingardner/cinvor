% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=24');


scenarioName = 'Tuning width modulation only';
testParameter = 'kappa';
testParameterValues = [15 20];


% display which scenario we are running
disp(scenarioName);
m = setCinvorModel(e,testParameter,testParameterValues(1));
sample = 1:10:180;
voxelResponse = m.ws*m.neuralResponse
voxelResponse = voxelResponse(:,sample);
F = fft(voxelResponse,[],2); %takes Fourier transform of each row.
avgTransform = fftshift(mean(abs(F),1));
avgTransform = min(1,avgTransform); %makes plots nicer; 
clf;
plot(avgTransform);
%plot(80:102,avgTransform(80:102));
%axis([80 102 0 1.5]);