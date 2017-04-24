preferredOrientations = 0:1:179;
angles = 0:1:179;
k = length(preferredOrientations);
n = length(angles);
x = linspace(0,2*pi,181);
x = x(1:end-1);


%rank 5 approximation rectified sin^2
for exponent = 1:10
	u = max(0,sin(x)).^exponent;
	F = abs(fftshift(fft(u)));
	plot(F(88:94));
	hold on
end
xlabel('Fourier Component');
ylabel('Magnitude');

