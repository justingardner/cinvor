trainKappa = [30, 15, 10, 5, 2];
h=mlrSmartfig(sprintf('results'));
for j = 1:length(trainKappa)
	kappa = [0.5,1,2,4,6,10,15,20,25,30,40,50];
	%kappa = [0.5,1,6];
	FWHM_vect = zeros(1,length(kappa));
	for i = 1:length(kappa)
		FWHM_vect(i) = vonMisesFWHM(kappa(i));
	end
	%kappa = [0.5,1,2,4,6];
	FWHM = zeros(2,length(kappa));
	for i = 1:length(kappa)
		FWHM(:,i) = simuTuningChange(1,trainKappa(j),kappa(i));
	end
	subplot(1,length(trainKappa),j);
	plot(FWHM_vect,FWHM(1,:),'bo-');
	hold on;
	plot(FWHM_vect,FWHM(2,:),'ro-');
	plot(FWHM_vect,FWHM_vect,'g')
	xlabel('neural tuning width (FWHM)','FontSize',13);
	legend('trained CRF','tested CRF','identity line');
	ylabel('channel response function width (FWHM)','FontSize',13);
	title(['decoded CRFs trained on FWHM = ',num2str(vonMisesFWHM(trainKappa(j)))],'FontSize',13);
	set(gca,'fontsize',13)
end