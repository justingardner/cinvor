channelPrefs = [0:22.5:180];
rangeScaleFac = 2;
% figure 3 (model fit to data)
thisData = load('simuCinvorData2');
allChannelTuning = thisData.allChannelTuning;
std_mat = thisData.std_mat;
realData = load('allData');
realChannelTuning2 = reshape(realData.allChannelTuning(2,:,:),[2,9]);
errorChannelTuning = reshape(realData.errorChannelTuning(2,:,:),[2,9]);
scenario=1;
scenarioName='Amplitude modulation only';

clf;
h=mlrSmartfig(sprintf('simuCinvor%i',scenario));
set(h,'name',sprintf('Scenario: %s',scenarioName));

%allChannelFit(:,1)=allChannelFit(:,1)*m(1).rangeScaleFac;
for i=1:2
  subplot(2,1,i);
  xp= channelPrefs(1):1:channelPrefs(end)*rangeScaleFac;
  %yp=vonMises(allChannelFit(i,1:4),d2r(xp));
  xpOri=xp/rangeScaleFac;
  l(1) = plot(channelPrefs, 1/2*(allChannelTuning(i,:) + fliplr(allChannelTuning(i,:))),['k']); hold on;
  errorbar(channelPrefs, 1/2*(allChannelTuning(i,:) + fliplr(allChannelTuning(i,:))),1.96*std_mat(i,:),['ko']);
  %plot(xpOri,yp,[getcolor(i),'--']);
  axis([-10,190,0,0.35]);
  hold on;
  l(2) = plot(channelPrefs,realChannelTuning2(i,:),'or','MarkerFaceColor','r','MarkerEdgeColor',[1,1,1],'MarkerSize',10);
  errbar(channelPrefs,realChannelTuning2(i,:),errorChannelTuning(i,:),'-r');
  drawPublishAxis;
  if(i == 1)
    hTitle = title('Channel response functions low contrast');
  else
    hTitle = title('Channel response functions high contrast');
  end
  hXLabel = xlabel('Orientation (deg)');
  hYLabel = ylabel('Channel response function');
  hLegend = legend([l(1),l(2)],{'Model','Data'});
  set( gca                       , ...
    'FontName'   , 'Helvetica' );
	set([hLegend,hTitle, hXLabel, hYLabel], ...
	    'FontName'   , 'AvantGarde');
	set([hLegend,gca]             , ...
	    'FontSize'   , 11           );
	set([hXLabel, hYLabel]  , ...
	    'FontSize'   , 12          );
	set( hTitle                    , ...
	    'FontSize'   , 14          , ...
	    'FontWeight' , 'bold'      );

	set(gca, ...
	  'Box'         , 'off'     , ...
	  'TickDir'     , 'out'     , ...
	  'TickLength'  , [.02 .02] , ...
	  'XMinorTick'  , 'on'      , ...
	  'YMinorTick'  , 'on'      , ...
	  'YGrid'       , 'on'      , ...
	  'XColor'      , [.3 .3 .3], ...
	  'YColor'      , [.3 .3 .3], ...
	  'YTick'       , 0:0.1:1, ...
	  'LineWidth'   , 1         );
end


%change tuning figure
thisData = load('widthChangeData3');
lowKappaVals = thisData.lowKappaVals;
allChannelTuning = thisData.allChannelTuning;

h=mlrSmartfig(sprintf('simuCinvor%i',scenario));
set(h,'name',sprintf('Scenario: %s',scenarioName));

%allChannelFit(:,1)=allChannelFit(:,1)*m(1).rangeScaleFac;
for i=1:2
  subplot(2,1,i);
  xp= channelPrefs(1):1:channelPrefs(end)*rangeScaleFac;
  %yp=vonMises(allChannelFit(i,1:4),d2r(xp));
  xpOri=xp/rangeScaleFac;
  legendVect = [];
  legendLabel = {};
  for jj = 1:length(lowKappaVals)
    l(jj) = plot(channelPrefs, reshape((allChannelTuning(jj,i,:)),[1,9]),'Color',[0.2,jj/length(lowKappaVals),1 - jj/length(lowKappaVals)]); hold on;
    legendVect = [legendVect, l(jj)];
    legendLabel{jj} = ['neural tuning fwhm ', num2str(vonMisesFWHM(lowKappaVals(jj)))];
  end
  %plot(xpOri,yp,[getcolor(i),'--']);
  axis([-10,190,0,0.40]);
  if(i == 1)
    hTitle = title('Channel response functions low contrast');
  else
    hTitle = title('Channel response functions high contrast');
  end
  hXLabel = xlabel('Orientation (deg)');
  hYLabel = ylabel('Channel response function');
  hLegend = legend(legendVect,legendLabel);
  set( gca                       , ...
    'FontName'   , 'Helvetica' );
  set([hLegend,hTitle, hXLabel, hYLabel], ...
      'FontName'   , 'AvantGarde');
  set([hLegend,gca]             , ...
      'FontSize'   , 11           );
  set([hXLabel, hYLabel]  , ...
      'FontSize'   , 12          );
  set( hTitle                    , ...
      'FontSize'   , 14          , ...
      'FontWeight' , 'bold'      );

  set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:0.1:1, ...
    'LineWidth'   , 1         );
end



%changeWidthFigure (sine channel)
thisData = load('widthSearch6');
kappaVals = thisData.kappaVals;
allWidths = thisData.allWidths;
allR2 = thisData.allR2;
color_array = [105,217,255;
140,232,121;
255,226,145;
232,129,121;
172,133,255];
figure;
legend_array = {};
for i = 1:length(kappaVals)
  toPlot = allR2(i,:) > 0;
  plot(allR2(i,toPlot),allWidths(i,toPlot),'Color',color_array(i,:)/256);
  hold on;
  xlabel('r2');;
  ylabel('crf fwhm');
  title('Change Neural Tuning Width (Sine Channels)');
  legend_array{i} = ['neural tuning width = ',num2str(vonMisesFWHM(kappaVals(i)))];

end

legend(legend_array);
axis([0,1,20,60]);
drawPublishAxis;

%changeWidthFigure (stick channel)
thisData = load('widthSearch5');
kappaVals = thisData.kappaVals;
allWidths = thisData.allWidths;
allR2 = thisData.allR2;
color_array = [105,217,255;
140,232,121;
255,226,145;
232,129,121;
172,133,255];
figure;
legend_array = {};
for i = 1:length(kappaVals)
  toPlot = allR2(i,:) > 0;
  plot(allR2(i,toPlot),allWidths(i,toPlot),'Color',color_array(i,:)/256);
  hold on;
  xlabel('r2');
  ylabel('crf fwhm')
  title('Change Neural Tuning Width (Stick Channels)');
  legend_array{i} = ['neural tuning width = ',num2str(vonMisesFWHM(kappaVals(i)))];
end

legend(legend_array);
axis([0,1,10,60]);
drawPublishAxis;



%changeChannelWidthFigure
thisData = load('channelWidthSearch');
expVals = thisData.expVals;
allWidths = thisData.allWidths;
allR2 = thisData.allR2;
color_array = [105,217,255;
140,232,121;
255,226,145;
232,129,121;
172,133,255];
figure;
legend_array = {};
for i = 1:length(expVals)
  toPlot = allR2(i,:) > 0;
  plot(allR2(i,toPlot),allWidths(i,toPlot),'Color',color_array(i,:)/256);
  hold on;
  xlabel('r2');
  ylabel('crf fwhm')
  title('Change Basis Width');
  legend_array{i} = ['channel basis exp = ',num2str(expVals(i))];
end

legend(legend_array);
axis([0,1,10,60]);
drawPublishAxis;

