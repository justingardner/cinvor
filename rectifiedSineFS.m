preferredOrientations = 0:1:179;
angles = 0:1:179;
k = length(preferredOrientations);
n = length(angles);
responses = zeros(k,n);
for i = 1:k
	for j = 1:n
		%responses(i,j) = 0.75*cos(pi*(angles(j) - preferredOrientations(i))/90)+0.25*cos(3*pi*(angles(j) - preferredOrientations(i))/90);
		responses(i,j) = cos(pi*(angles(j) - preferredOrientations(i))/90).^2;
		%myAngle = mod(angles(j) - preferredOrientations(i),180);
		%if(myAngle > 90)
		%	myAngle = myAngle - 180;
		%end
		%responses(i,j) = exp(-(pi*(myAngle)/90)^2);
	end
end

rangel = 60; %plotting range
rangeu = 120;
h=mlrSmartfig(sprintf('plots'));

t = linspace(-pi,pi,180);
Fs = 2*pi/180; %sampling frequency
f = linspace(-90.5*Fs,89.5*Fs,180);
f = f(rangel:rangeu);
subplot(4,2,1);
plot(t,responses(90,:))
F = fft(responses(90,:))
xlabel('t');
ylabel('cos(t)');
title('Cosine Squared');
axis([-pi pi -1 1])
subplot(4,2,2);
Fplot = fftshift(abs(F));
plot(f,Fplot(rangel:rangeu));
xlabel('s');
ylabel('F(s)');
title('Fourier Transform');
axis([-1.5 1.5 -100 100]);
responses(90,[1:45, 136:180]) = 0;
subplot(4,2,5);
plot(t,responses(90,:))
xlabel('t');
ylabel('max(0,cos(t))');
title('Rectified Cosine Squared');
axis([-pi pi -1 1])
F = fft(responses(90,:))
subplot(4,2,6);
Fplot = fftshift(abs(F));
plot(f,Fplot(rangel:rangeu));
xlabel('s');
ylabel('F(s)');
title('Fourier Transform');
axis([-1.5 1.5 -100 100]);
F(logical(abs(fft(responses(90,:))) < 3)) = 0;
subplot(4,2,8);
Fplot = fftshift(abs(F));
plot(f,Fplot(rangel:rangeu));
xlabel('s');
ylabel('F(s)');
title('Fourier Transform');
axis([-1.5 1.5 -100 100]);
subplot(4,2,7);
plot(t,ifft(F));
xlabel('t');
ylabel('f(t)');
title('Rank 7 Steerable Approximation of Rectified Cosine Squared');
axis([-pi pi -1 1])
newFunc = ifft(F);
for i = 1:180
	M(i,:) = circshift(newFunc',i)';
end
rank(M); %7

subplot(4,2,3);
plot(t,responses(90,:)>0)
F = fft(responses(90,:)>0)
xlabel('t');
ylabel('box(t)');
title('Box Function');
axis([-pi pi -1 1])
subplot(4,2,4);
Fplot = fftshift(abs(F));
plot(f,Fplot(rangel:rangeu));
xlabel('s');
ylabel('F(s)');
title('Fourier Transform');
axis([-1.5 1.5 -100 100]);
