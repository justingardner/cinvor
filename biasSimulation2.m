%Notation follows The Prediction Properties of Inverse and Reverse Regression 
%for the Simple Linear Calibration Problem

total = zeros(1,21);
T = 10000;

for i = 1:T
	B_0 = 0;
	B_1 = 1;
	sigma = 2;

	n = 100;
	x = zeros(n,1);
	x(n/2+1:end) = 1;
	epsilon = randn(n,1)*sigma;
	y = B_0 + B_1*x + epsilon;
	y_mean = mean(y);
	y_star = y - y_mean;
	x_mean = mean(x);

	B_1_hat = sum((x-x_mean).*y_star)/sum((x-x_mean).*(x - x_mean));
	B_0_hat_star = -B_1_hat*x_mean;
	B_0_hat = B_0_hat_star + y_mean;

	y_star_hat = B_0_hat_star + B_1_hat*x;
	sigma_sq_star = sum((y_star - y_star_hat).^2)/(n-2);

	x_to_pred = -10:10;
	y_0 = B_1*x_to_pred + B_0;
	x_0 = (y_0 - B_0_hat)/B_1_hat;

	bias_factor = sigma^2/(sum((x-x_mean).*(x - x_mean)));
	total = total + 1/T*x_0;
end


%bias away from mean!
plot(x_to_pred,total,'o');
hold on;
plot(x_to_pred,x_to_pred);