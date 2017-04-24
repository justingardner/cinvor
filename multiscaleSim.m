x = linspace(0,2*pi,180)
u = max(0,6 - 0.25*(x - pi).^4);
F = fft(u);
F(logical(abs(real(F)) < 20 | abs(real(F)) > 3000)) = 0;
%plot(real(fftshift(F)));
u = ifft(F);
plot(u);
U = zeros(6,180);
for i = 1:180
	U(1,i) = u(i) + u(mod(-i - 1,180) + 1);
	U(2,i) = u(i) - u(mod(-i - 1,180) + 1);
	U(3,i) = u(mod(i - 60 - 1,180)+1) + u(mod(-i - 60 - 1,180) + 1);
	U(4,i) = u(mod(i - 60 - 1,180)+1) - u(mod(-i - 60 - 1,180) + 1);
	U(5,i) = u(mod(i - 120 - 1,180)+1) + u(mod(-i - 120 - 1,180) + 1);
	U(6,i) = u(mod(i - 120 - 1,180)+1) - u(mod(-i - 120 - 1,180) + 1);
end
p=@(x) sign(x).*log(abs(x).*(exp(pi)-1)/pi+1);

V = zeros(5,180);
for i = 1:180
	for n = 1:5
		V(n,i) = interp1(x,u,mod(p(x(i)) -n*pi/7,2*pi),'spline') + interp1(x,u,mod(-p(x(i)) -n*pi/7,2*pi),'spline');
	end
end
newRow = [fliplr(V(5,:)) V(5,:)];
shift = 50;
V = [V; newRow(181-shift:360-shift)];