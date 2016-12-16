
%noise,kappa,alpha
ROI = 2;
M = 8;
data = load('channelTuning');
realChannelTuning = reshape(data.totalChannelTuning(ROI,:,:),[2,M+1]);
f = @(x)simuCinvorPlot(x(1),x(2),x(3),realChannelTuning);
x_base = [8,2.5,0.5];

param_vals = [3,7,11;
			  1,2.5,5;
			  0.2,0.5,1];%[1,3,5,7,9,11,13];
%x = fmincon(f,x_0,A,b);
h=mlrSmartfig(sprintf('simuCinvor%i',scenario));
set(h,'name',sprintf('Scenario: %s',scenarioName));
channelPrefs = [e.stimVals 180];
param_names = {'noise','kappa','alpha'}

allChannelTuning = zeros(length(x_base),size(param_vals,1),length(testParameterValues),M+1);
for pNumber = 1:length(x_base)
	
	for j = 1:size(param_vals,1)
		x = x_base;
		x(pNumber) = param_vals(pNumber,j);
		allChannelTuning(pNumber,j,:,:) = reshape(f(x),size(allChannelTuning(pNumber,j,:,:))); 
	end


	%allChannelFit(:,1)=allChannelFit(:,1)*m(1).rangeScaleFac;
	for i=1:length(testParameterValues)
	  subplot(2,3,3*(i-1)+pNumber);
	  xp= channelPrefs(1):1:channelPrefs(end)*m(1).rangeScaleFac;
	  %yp=vonMises(allChannelFit(i,1:4),d2r(xp));
	  xpOri=xp/m(1).rangeScaleFac;
	  for jj = 1:length(noise_vals)
	    plot(channelPrefs, reshape((allChannelTuning(pNumber,jj,i,:)),[1,9]),getcolor(jj)); hold on;
	  end
	  %plot(xpOri,yp,[getcolor(i),'--']);
	  axis([-10,190,-0.05,0.25]);
	  hold on;
	  plot(channelPrefs,realChannelTuning(i,:),'g');
	  if(i == 1)
	    title([param_names{pNumber},' : Channel tuning functions low contrast']);
	  else
	    title([param_names{pNumber},' : Channel tuning functions high contrast']);
	  end
	end
end