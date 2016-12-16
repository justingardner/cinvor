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
u = max(0,sin(x)).^2;
F = fft(u);
F(logical(abs(F) < 10)) = 0;
u = ifft(F)

%u = sin(x).^2;
u = sin(x/2);

%rectified sin^2
u = max(0,sin(x)).^3;

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
h=mlrSmartfig(sprintf('plots'));
mycolor = rand(8,3);
subplot(4,1,1);
I = [1:23:180]; %[10,30,50,70,90,110];
j = 1;
for index = I
	plot(responses(index,:),'color',mycolor(j,:));
	hold on;
	j = j + 1;
end
xlabel('angle (degrees)');
ylabel('function value');
title('rectified sin^2 basis functions');

subplot(4,1,2);
index = 10;
I = [1:23:180]; %[10,30,50,70,90,110];
M = responses(I,:);
thisVect = responses(index,:) + responses(index+90,:);
%thisVect = responses(index,:);
w = M'\thisVect';


plot(M'*w,'r')
hold on
plot(thisVect,'b')
legend('best fit from basis','original shifted sine squared');
xlabel('angle (degrees)');
ylabel('function value');
title('shifted sin^2 can be reached by this basis');

subplot(4,1,3);
index = 10;
I = [1:23:180]; %[10,30,50,70,90,110];
M = responses(I,:);
thisVect = responses(index,:);
%thisVect = responses(index,:);
w = M'\thisVect';


plot(M'*w,'r')
hold on
plot(thisVect,'b')
legend('best fit from basis','original shifted basis');
xlabel('angle (degrees)');
ylabel('function value');
title('shifted, rectified sin^2 cannot be reached by this basis');

subplot(4,1,4);

I = [1:15:180]; %[10,30,50,70,90,110];
j = 1;
for index = I
	plot(w(j)*responses(index,:),'color',mycolor(j,:));
	hold on;
	j = j+1;
end
xlabel('angle (degrees)');
ylabel('function value');
title('best fit of basis functions to shifted, rectified sin^2');

%rank 5 approximation rectified sin^5
u = max(0,sin(x)).^5;
F = fft(u);
F(logical(abs(F) < 15)) = 0;
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

I = [1:15:180]; %[10,30,50,70,90,110];
M = responses(I,:);
index = 31;
w = M'\responses(index,:)';

%subplot(2,1,2);
%plot(M'*w,'r')
%hold on
%plot(responses(index,:),'b')
%xlabel('angle (degrees)');
%ylabel('function value');
%title('rank 5 steerable approximation of sin^5');