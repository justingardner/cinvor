responses = reshape(testChannelResponse,192*8,1);
channelVals = [0,0.0883883476483185,1];
threshold = 0.01;
predictedVals = zeros(1,length(channelVals));
for i = 1:length(channelVals)
	predictedVals(1,i) = mean(responses(logical(abs(reshape(channel.channelResponse,192*8,1) - channelVals(i)) < threshold)));
end
figure;
plot(channelVals,predictedVals);
hold;
plot(channelVals,channelVals,'r');
xlabel('actual channel acitvation');
ylabel('predicted channel acitvation');
title('actual vs predicted channel activation');

figure;
plot(channel.idealStimVals,avgChannelResponse);
hold;
plot(channel.idealStimVals,channel.idealStimResponse(5,:),'r');