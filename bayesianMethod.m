noise = 0.2;
total = zeros(2,2);
avg_pred = 0;
T = 1000;
n_voxels = 2;
for j = 1:T
	x = [1, 0; 0, 1];
	true_w = rand(n_voxels,2);
	true_w = [1 0; 0 1];
	true_w(1,1) = rand()*0.3+0.7;
	true_w(2,2) = true_w(1,1);
	true_w(2,1) = 1 - true_w(1,1);
	true_w(1,2) = 1 - true_w(1,1);
	true_y = true_w*x;
	y_obs = randn(n_voxels,2)*noise+ true_w*x;
	w_hat = y_obs*x'*inv(x*x');
	x_hat = inv(w_hat'*w_hat)*w_hat'*true_y;
	total = total + 1/T*x_hat;
	hypothesis = linspace(0.7,1,100);
	likelihood = zeros(1,length(hypothesis));
	value = zeros(1,length(hypothesis));
	for i = 1:length(hypothesis)
		h = hypothesis(i);
		my_w = [h, 1-h; 1-h, h];
		pred_y = my_w*x;
		likelihood(i) = exp(-sum(sum((pred_y - y_obs).^2))/(2*noise^2));
		channel_pred = inv(my_w)*true_y;
		value(i) = channel_pred(1,1);
	end
	likelihood = likelihood/sum(likelihood);
	a = sum(likelihood.*value);
	avg_pred = avg_pred + 1/T*sum(likelihood.*value);
end


x = [1, 0; 0, 1];
true_w = rand(n_voxels,2);
true_w = [1 0; 0 1];
true_w(1,1) = rand()*0.3+0.7;
true_w(2,2) = true_w(1,1);
true_w(2,1) = 1 - true_w(1,1);
true_w(1,2) = 1 - true_w(1,1);
hypothesis = linspace(0.7,1,100);
likelihood = zeros(1,length(hypothesis));
value = zeros(1,length(hypothesis));
for i = 1:length(hypothesis)
	h = hypothesis(i);
	my_w = [h, 1-h; 1-h, h];
	pred_y = my_w*x;
	likelihood(i) = exp(-sum(sum((pred_y - y_obs).^2))/(2*noise^2));
	channel_pred = inv(my_w)*true_y;
	value(i) = channel_pred(1,1);
end
likelihood = likelihood/sum(likelihood)
avg_pred = sum(likelihood.*value);


%doesn't work
total = zeros(1,1);
T = 10000;
trial = zeros(1,T);
for j = 1:T
	true_w = rand()*0.3+0.7;
	x = [1];
	true_y = true_w*x;
	y_obs = randn(1,1)*noise+ true_w*x;
	w_hat = y_obs*x'*inv(x*x');
	x_hat = inv(w_hat'*w_hat)*w_hat'*true_y;
	hypothesis = linspace(0.7,1,100);
	likelihood = zeros(1,length(hypothesis));
	value = zeros(1,length(hypothesis));
	for i = 1:length(hypothesis)
		h = hypothesis(i);
		my_w = [h];
		pred_y = my_w*x;
		likelihood(i) = exp(-(pred_y - y_obs).^2/(2*noise^2));
		channel_pred = inv(my_w)*true_y;
		value(i) = channel_pred;
	end
	likelihood = likelihood/sum(likelihood);
	avg_pred = sum(likelihood.*value);
	total = total + 1/T*avg_pred;
	trial(j) = avg_pred;
end


total = zeros(1,1);
T = 10000;
trial = zeros(1,T);
for j = 1:T
	w_param = rand()*0.3+0.7;
	true_w = [w_param, (1-w_param); (1-w_param), w_param];
	x = [1, 0; 0, 1];
	true_y = true_w*x;
	y_obs = randn(2,2)*noise+ true_w*x;
	w_hat = y_obs*x'*inv(x*x');
	x_hat = inv(w_hat'*w_hat)*w_hat'*true_y;
	hypothesis = linspace(0.7,1,100);
	likelihood = zeros(1,length(hypothesis));
	value = zeros(1,length(hypothesis));
	for i = 1:length(hypothesis)
		h = hypothesis(i);
		my_w = [h];
		pred_y = my_w*x;
		likelihood(i) = exp(-(pred_y(1,1) - y_obs(1,1)).^2/(2*noise^2))*exp(-(pred_y(2,1) - y_obs(2,1)).^2/(2*noise^2))*exp(-(pred_y(2,2) - y_obs(2,2)).^2/(2*noise^2))*exp(-(pred_y(1,2) - y_obs(1,2)).^2/(2*noise^2));
		channel_pred = inv(my_w)*true_y;
		value(i) = channel_pred(1,1);
	end
	likelihood = likelihood/sum(likelihood);
	avg_pred = sum(likelihood.*value);
	total = total + 1/T*avg_pred;
	trial(j) = avg_pred;
end