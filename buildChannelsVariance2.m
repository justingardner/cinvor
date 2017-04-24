global myRealChannel taoParam sigmaParam rhoParam;
myRealChannel = realChannel; %channel fit from true data (only channelWeights will be used)
taoParam = 0.08;
sigmaParam = 0.02;
rhoParam = 0.15;
N_cond = 1;
rho_mean = zeros(N_cond,1);
sigma_mean = zeros(N_cond,1);
tao_mean = zeros(N_cond,length(realChannel.tao));
rho_var = zeros(N_cond,1);
sigma_var = zeros(N_cond,1);
tao_var = zeros(N_cond,length(realChannel.tao));
%rhoCond = [0.15,.22 , .05, .1, .2];
%sigmaCond = [0.02, .59, .03, .7, .2];
%taoCond = [0.08, .82, .06, .3, .9];
rhoCond = [.16];
sigmaCond = [1.02];
taoCond = [mean(myRealChannel.tao)];
% set the experiment (number of orientations to test and repetitions)
e = specifyCinvorExperiment('stimLevel=8','trialPerStim=24');
m = setCinvorModel(e);
trainInstances = getCinvorInstances(m,e); %these don't matter, just to pass in something.
T = 3;
for condition = 1:N_cond
	rhoParam = rhoCond(condition);
	%taoParam = taoCond(condition);
	taoParam = myRealChannel.tao;
	sigmaParam = sigmaCond(condition);
	rho_list = zeros(T,1);
	sigma_list = zeros(T,1);
	tao_list = zeros(T,length(realChannel.tao));
	for timeTrial = 1:T
		channel = buildChannelsTesting(trainInstances,e.stimVals);
		rho_list(timeTrial) = channel.rho;
		sigma_list(timeTrial) = channel.sigma;
		tao_list(timeTrial,:) = channel.tao;
	end
	rho_mean(condition) = mean(rho_list);
	sigma_mean(condition) = mean(sigma_list);
	tao_mean(condition,:) = mean(tao_list,1);
	rho_var(condition) = var(rho_list);
	sigma_var(condition) = var(sigma_list);
	tao_var(condition,:) = var(tao_list,1);
end

mean_tao_mean = mean(tao_mean,2);
var_tao_mean = 1/size(tao_var,2)*mean(tao_var,2);


figure;
subplot(1,3,1);
plot(1:N_cond,rhoCond,'bo');
hold on;
plot(1:N_cond,rho_mean,'ro')
errorbar(1:N_cond,rho_mean,2*sqrt(rho_var),'r');
title('rho');

subplot(1,3,2);
plot(1:N_cond,sigmaCond,'bo');
hold on;
plot(1:N_cond,sigma_mean,'ro')
errorbar(1:N_cond,sigma_mean,2*sqrt(sigma_var),'r');
title('sigma');

subplot(1,3,3);
plot(1:N_cond,taoCond,'bo');
hold on;
plot(1:N_cond,mean_tao_mean,'ro')
errorbar(1:N_cond,mean_tao_mean,2*sqrt(var_tao_mean),'r');
title('tao');