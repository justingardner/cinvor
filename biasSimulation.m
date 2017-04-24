%change noise
noise_val = 0:0.01:1;
low_val = zeros(1,length(noise_val));
high_val = zeros(1,length(noise_val));
for j = 1:length(noise_val)
	noise = noise_val(j);
	total = zeros(2,2);
	T = 1000;
	n_voxels = 100;
	for i = 1:T
		x = [100, 99; 99, 100];
		true_w = rand(n_voxels,2);
		true_y = true_w*x;
		y_obs = randn(n_voxels,2)*noise+ true_w*x;
		w_hat = y_obs*x'*inv(x*x');
		x_hat = inv(w_hat'*w_hat)*w_hat'*true_y;
		total = total + 1/T*x_hat;
	end
	low_val(j) = 1/2*(total(1,2) + total(2,1));
	high_val(j) = 1/2*(total(1,1) + total(2,2));
end
clf;
plot(noise_val,low_val,'o');
hold on;
plot(noise_val,high_val,'ro');
legend('estimate of c_{1,1}','estimate of c_{2,2}');
xlabel('standard deviation of noise','FontSize', 16);
ylabel('percent of total channel activation','FontSize', 16);
set(gca,'fontsize',14);

%change n_voxels
voxel_val = 2:100;
low_val = zeros(1,length(voxel_val));
high_val = zeros(1,length(voxel_val));
for j = 1:length(voxel_val)
	n_voxels = voxel_val(j);
	total = zeros(2,2);
	T = 1000;
	noise = 0.3;
	for i = 1:T
		x = [1, 0; 0, 1];
		true_w = rand(n_voxels,2);
		true_y = true_w*x;
		y_obs = randn(n_voxels,2)*noise+ true_w*x;
		w_hat = y_obs*x'*inv(x*x');
		x_hat = inv(w_hat'*w_hat)*w_hat'*true_y;
		total = total + 1/T*x_hat;
	end
	low_val(j) = 1/2*(total(1,2) + total(2,1));
	high_val(j) = 1/2*(total(1,1) + total(2,2));
end
figure;
plot(voxel_val,low_val,'o');
hold on;
plot(voxel_val,high_val,'ro');

%no bias here:
plot(true_w(:),w_hat(:),'o');
true_w(:)\w_hat(:); %about 1

%still, no bias here
big_second = true_w(:,2) > 0.9;
plot(true_w(big_second,1),w_hat(big_second,1));
true_w(big_second,1)\w_hat(big_second,1); %about 1





%regularization
total = zeros(2,2);
T = 10000;
n_voxels = 100;
lambda = 0.5;
for i = 1:T
	x = [1, 0; 0, 1];
	true_w = rand(n_voxels,2);
	true_y = true_w*x;
	noise = 0.3;
	y_obs = randn(n_voxels,2)*noise+ true_w*x;
	w_hat = y_obs*x'*inv(x*x' + lambda*eye(2));
	x_hat = inv(w_hat'*w_hat)*w_hat'*true_y;
	total = total + 1/T*x_hat;
end




%reverse regression :0  also seems to exhibit bias.
total = zeros(2,2);
T = 10000;
for i = 1:T
	n_voxels = 20;
	true_w = rand(2,n_voxels);
	x = [1, 0; 0, 1];
	true_y = pinv(true_w'*true_w)*true_w'*x;
	noise = 0.5;
	y_obs = randn(n_voxels,2)*noise+ true_y;
	w_hat = x*y_obs'*pinv(y_obs*y_obs');
	x_hat = w_hat*true_y;
	total = total + 1/T*x_hat;
end


%Simultaneously fit forward and reverse model
noise = 0.3;
total = zeros(2,2);
T = 1000;
n_voxels = 100;
x = [1, 0; 0, 1];
true_w = rand(n_voxels,2);
true_y = true_w*x;
y_obs = randn(n_voxels,2)*noise+ true_w*x;
w_hat = y_obs*x'*inv(x*x');
x_hat = inv(w_hat'*w_hat)*w_hat'*true_y;
v_hat = x*pinv(y_obs);
x_pred = v_hat*y_obs;
%it holds that v_hat*w_hat = ident, but not w_hat*v_hat is way off.
mean(mean(abs(w_hat*v_hat*true_y - true_y)));


