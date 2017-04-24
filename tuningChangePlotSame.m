trainKappa = [30, 15, 10, 5, 2];
neuralWidth = zeros(length(trainKappa),1);
for i = 1:length(trainKappa)
	neuralWidth(i) = vonMisesFWHM(trainKappa(i));
end
noise = [0,0.1,0.3,0.5,1,2,5];
FWHM = zeros(2,length(trainKappa),length(noise));
trials = 10;
for j = 1:length(trainKappa)
	for noiseIndex = 1:length(noise)
		for trial = 1:trials
			thisNoise = noise(noiseIndex);
			FWHM(:,j,noiseIndex) = FWHM(:,j,noiseIndex) + 1.0/trials*simuTuningChange(1,trainKappa(j),trainKappa(j),thisNoise);
		end
	end
end
clf;
for j = 1:length(trainKappa)
	plot(noise,reshape(FWHM(2,j,:),[length(noise),1]),getcolor(j));
	hold on;
end
xlabel('noise (noise to signal ratio)');
ylabel('channel function width FWHM)');
legend(cellstr(num2str(neuralWidth, 'neural width=%-d')))
