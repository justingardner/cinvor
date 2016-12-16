preferredOrientations = 0:1:179;
angles = 0:1:179;
k = length(preferredOrientations);
n = length(angles);
x = linspace(0,2*pi,181);
x = x(1:end-1);

%rectified sin^5
u = max(0,sin(x)).^5;

%rank 5 approximation rectified sin^5
u = max(0,sin(x)).^5;
F = fft(u);
F(logical(abs(F) < 15)) = 0;
u = ifft(F)

%rectified sin^2
u = max(0,sin(x)).^2;

%rank 5 approximation rectified sin^2
u = max(0,sin(x)).^3;
F = fft(u);
sortedComponents = sort(abs(F))
numToKeep = 7;
F(logical(abs(F) < sortedComponents(end - numToKeep+1))) = 0;
u = ifft(F)

%rank 5 approximation rectified sin^2
u = max(0,sin(x)).^5;
F = abs(fftshift(fft(u)));
sortedComponents = sort(abs(F));
numToKeep = 7;

F(logical(abs(F) < sortedComponents(end - numToKeep+1))) = 0;
u = ifft(F)


responses = zeros(k,n);
for i = 1:k
	for j = 1:n
		%responses(i,j) = 0.75*cos(pi*(angles(j) - preferredOrientations(i))/90)+0.25*cos(3*pi*(angles(j) - preferredOrientations(i))/90);
		%responses(i,j) = max(0,cos(pi*(angles(j) - preferredOrientations(i))/90))^5;
		%responses(i,j) = sin(p(x(i))-x(j));
		responses(i,j) = u(mod(i+j,180)+1);
		%myAngle = mod(angles(j) - preferredOrientations(i),180);
		%if(myAngle > 90)
		%	myAngle = myAngle - 180;
		%end
		%responses(i,j) = exp(-(pi*(myAngle)/90)^2);
	end
end

I = [1:23:180]; %[10,30,50,70,90,110];
M = responses(I,:);
index = 100;
w = M'\responses(index,:)';
