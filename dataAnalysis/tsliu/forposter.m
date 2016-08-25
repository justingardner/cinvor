

%figure 1
figure;
subplot(1,2,1);
x = linspace(0,2*pi,360);
plot(max(-cos(x),0).^2+0.05,'r','LineWidth',4);

subplot(1,2,2);
x = linspace(0,2*pi,360);
plot(max(-cos(x),0).^10+0.05,'r','LineWidth',4);

%figure 2
figure;
colorList = {'r','y',[0.25,0.8,0.25],[0.5,1,0.5],'c','b',[0.5,0,0.5],'m'};
for i = 1:8
	subplot(8,1,i);
	plot(channel.spanResponse(:,i)+0.05,'Color',colorList{i});
	axis([0 180 0 1.1]);
end

%figure 3
myMatrix = [1 0 1 0 0 0 0 0;
			0 0 0 0.5 0.5 0.5 0 0;
			0.7 0.7 0 0 0 0 1 0 ;
			0 0 0 1 0 0 1 0;
			0.4 0.5 0.6 0.6 0.5 0.6 0.4 0.5;
			1 0 1 0 1 0 1 0;
			1 0 0 0 0 0 0 0;
			0.7 0 0.7 0.7 0.7 0 0 0];

voxelResponse = myMatrix*channel.spanResponse';
figure;
for i = 1:8
	subplot(8,1,i);
	plot(voxelResponse(i,:)+0.05,'k');
	axis([0 180 0 1.1]);
end

%figure 4;
figure;
x = linspace(0,2*pi,360);
plot(max(-cos(x),0).^2+0.05,'k','LineWidth',1.5); hold on;
plot(0.5*max(-cos(x),0).^2+0.05,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
legend('high contrast','low contrast');
axis([0 360 0 1.1]);

%figure 5
figure;
colorList = {'r','y',[0.25,0.8,0.25],[0.5,1,0.5],'c','b',[0.5,0,0.5],'m'};
for i = 1:8
	subplot(8,1,i);
	plot(channel.spanResponse(:,i)+0.05,'Color',colorList{i}); hold on;
	plot(0.5*channel.spanResponse(:,i)+0.05,':','Color',colorList{i});
	axis([0 180 0 1.1]);
end

%figure 5
myMatrix = [1 0 1 0 0 0 0 0;
			0 0 0 0.5 0.5 0.5 0 0;
			0.7 0.7 0 0 0 0 1 0 ;
			0 0 0 1 0 0 1 0;
			0.4 0.5 0.6 0.6 0.5 0.6 0.4 0.5;
			1 0 1 0 1 0 1 0;
			1 0 0 0 0 0 0 0;
			0.7 0 0.7 0.7 0.7 0 0 0];

voxelResponse = myMatrix*channel.spanResponse';
figure;
for i = 1:8
	subplot(8,1,i);
	plot(voxelResponse(i,:)+0.05,'k'); hold on;
	plot(0.5*voxelResponse(i,:)+0.05,':','Color','k');
	axis([0 180 0 1.1]);
end


%figure 6
figure;
colorList = {'r','y',[0.25,0.8,0.25],[0.5,1,0.5],'c','b',[0.5,0,0.5],'m'};
for i = 1:8
	subplot(8,1,i);
	plot(channel.spanResponse(:,i)+0.05,'Color',colorList{i}); hold on;
	plot(channel.spanResponse(:,i).^5+0.05,':','Color',colorList{i});
	axis([0 180 0 1.1]);
end

%figure 6
myMatrix = [1 0 1 0 0 0 0 0;
			0 0 0 0.5 0.5 0.5 0 0;
			0.7 0.7 0 0 0 0 1 0 ;
			0 0 0 1 0 0 1 0;
			0.4 0.5 0.6 0.6 0.5 0.6 0.4 0.5;
			1 0 1 0 1 0 1 0;
			1 0 0 0 0 0 0 0;
			0.7 0 0.7 0.7 0.7 0 0 0];

voxelResponse = myMatrix*channel.spanResponse';
voxelResponse2 = myMatrix*(channel.spanResponse').^5;
figure;
for i = 1:8
	subplot(8,1,i);
	plot(voxelResponse(i,:)+0.05,'k'); hold on;
	plot(voxelResponse2(i,:)+0.05,':','Color','k');
	axis([0 180 0 1.1]);
end