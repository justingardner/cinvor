ROI = 2;
h=mlrSmartfig(sprintf('simuCinvor%i',scenario));
set(h,'name',sprintf('Scenario: %s',scenarioName));
subplot(4,2,1);
plot(e.stimVals,trainSTD(:,:,1),'o');
hold;
plot(e.stimVals,realStats.trainSTD(ROI,:,1),'ro');
title('trial to trial std (low contrast)');
legend('model','data');
subplot(4,2,2);
plot(e.stimVals,trainSTD(:,:,2),'o');
hold;
plot(e.stimVals,realStats.trainSTD(ROI,:,2),'ro');
title('trial to trial std (low contrast)');
legend('model','data');
subplot(4,2,3);
plot(e.stimVals,trainTestCorr(:,:,1),'o');
hold;
plot(e.stimVals,realStats.trainTestCorr(ROI,:,1),'ro');
title('correlation of train and test instances (low contrast)');
legend('model','data');
subplot(4,2,4);
plot(e.stimVals,trainTestCorr(:,:,2),'o');
hold;
plot(e.stimVals,realStats.trainTestCorr(ROI,:,2),'ro');
title('correlation of train and test instances (high contrast)');
legend('model','data');

subplot(4,2,5:6)
plot(channelWeightMean,'o');
hold;
plot(realStats.channelWeightMean(ROI,:),'ro');
title('mean of channel weights');
legend('model','data');

subplot(4,2,7:8)
plot(channelWeightStd,'o');
hold;
plot(realStats.channelWeightStd(ROI,:),'ro');
title('std of channel weights');
legend('model','data');