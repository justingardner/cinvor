M = 16;
avg = zeros(1,M);
avgVar = zeros(1,M);
T = 10;
recordings = zeros(T,M);
for ii = 1:T
	simuVoxelFourierTrialFN;
	recordings(ii,:) = plotVect;
	avg = avg + plotVect;
	%avgVar = avgVar + totalVar;
end
avg = avg/T;
%avgVar = avgVar/(T^2);
avgSTD = sqrt(avgVar);
clf
plot(avg);
hold on;
plot(groundTruth,'r');
%errorbar(1:M,avg,avgSTD*1.96)
%axis([0 M+1 0 5]);

%verify correct calculation of errorbars
%avgVar = avgVar*T; 
recordingsVar = N/(N-1)*moment(recordings,2,1);