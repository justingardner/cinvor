%noise,kappa,alpha
ROI = 2;
M = 8;
data = load('channelTuning');
realChannelTuning = reshape(data.totalChannelTuning(ROI,:,:),[2,M+1]);
f = @(x)simuCinvorSearch(x(1),x(2),x(3),realChannelTuning);
x_vals = [8,2.5,0.5;
       4, 5, 0.5;
       8, 10, 0.2;
       15,5, 0.2;
       15,5, 1;
       10,2,1];
A = [-1 0 0; 0 -1 0; 0 0 -1; 0 0 1];
b = [0; 0; 0; 1];
x_results = zeros(6,3);
%x = fmincon(f,x_0,A,b);
for iter = 1:6
	x_0 = x_vals(iter,:)
	x = patternsearch(f,x_0,A,b)
	x_results(iter,:) = x;
end

