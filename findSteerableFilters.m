preferredOrientations = 0:1:179;
angles = 0:1:179;
k = length(preferredOrientations);
n = length(angles);
x = linspace(0,2*pi,181)
x = x(1:end-1);
p=@(x) sign(x).*log(abs(x).*(exp(pi)-1)/pi+1);
u = sin(p(x));
responses = zeros(k,n);
for i = 1:k
	for j = 1:n
		%responses(i,j) = 0.75*cos(pi*(angles(j) - preferredOrientations(i))/90)+0.25*cos(3*pi*(angles(j) - preferredOrientations(i))/90);
		%responses(i,j) = max(0,cos(pi*(angles(j) - preferredOrientations(i))/90))^2;
		%responses(i,j) = sin(p(x(i))-x(j));
		myAngle = mod(angles(j) - preferredOrientations(i),180);
		%if(myAngle > 90)
		%	myAngle = myAngle - 180;
		%end
		%responses(i,j) = exp(-(pi*(myAngle)/90)^2);
	end
end