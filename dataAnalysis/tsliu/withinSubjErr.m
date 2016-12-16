% Program for sizing error bars when dealing with witin-subjects 
% experimental designs. Assumes data is set up with each participant 
% as a row, and each column as a condition.  
% Written by James Miller modified by Taosheng Liu  Michigan State University 2013
% Gozde S. % "I got 99 problems but within-subjects confidence intervals aint one."
% Reference: Cousineau 2005 "Confidence intervals in within-subject design: 
% A simpler solution to Loftus and Masson's method" 
% Tutorials in Quantitative Methods for Psychology, 2005, vol 1, 42-45.

function newErrors = withinSubjErr(data)
disp('Calculating within-subject error bars...');
disp('Assuming each subject is in a row and each column is a condition');
nsubj = size(data, 1); ncond = size(data, 2); newYs = zeros(nsubj, ncond);
subjMeans = mean(data, 2); condMeans = mean(data);
for i = 1:ncond
  for j = 1:nsubj
    newYs(j,i) = data(j,i)-subjMeans(j)+condMeans(i);
  end
end

newErrors = std(newYs)/sqrt(nsubj);
return;


%Test code
% Below you'll find the example data that was in the Cousineau paper cited above.  Evaluate the code in this section to see the result.
% The first figure that plots should match page 43 and the second figure should
% match theirs on page 45 except in that the first plot has the x values for bottom 
% line is barely shifted purely for demonstrative purposes so the errorbars don't overlap.  
% Location of the data file:
% http://www.tqmp.org/Content/vol01-1/p042/p042.dat

  data1 = [496.6435  496.9672  496.8114  489.6242  484.6379; 513.0236  538.2246  534.0575  512.0894  526.4623; 547.3417  542.4835  536.9625  528.6456  516.2566;...
  565.9576  554.6250  563.5193  542.8214  553.5758; 580.3726  574.6816  590.1554  586.2111  578.6230; 606.6820  606.1392  599.7016  595.9525  619.2350;...
  604.2611  626.7651  641.6801  640.7080  618.3975; 633.1669  639.0650  650.5923  634.4741  645.8410; 665.8995  664.1079  663.6813  666.0743  643.5422;...
  676.5453  688.5991  690.7779  674.4994  692.0214; 694.3974  712.3280  701.4573  689.8099  685.5879; 719.6566  727.1880  736.8660  724.8033  725.0685;...
  741.3328  742.7191  749.8870  740.8231  743.8731; 762.8086  765.7397  764.1394  760.9007  774.0160; 791.9821  792.8631  776.0112  787.8207  775.5864;...
  802.3982  818.2564  805.1859  800.6814  789.6591];

  data2 = [503.5429  477.8233  484.0521  476.2517  473.8266; 507.8191  496.8105  495.9390  506.8325  479.3559; 523.0459  543.4251  526.9505  517.1195  520.0745;...
  552.1515  563.9514  565.3441  543.5078  534.3435; 585.3215  583.3521  577.3410  584.6055  549.0680; 589.4510  575.8365  596.7620  577.2900  567.6478;...
  619.3001  621.1451  614.1312  597.6707  583.5788; 631.6709  618.1561  608.7894  614.4030  603.2316; 660.7510  644.0570  646.3204  627.4396  655.5683;...
  665.8861  678.1308  668.4901  659.3903  654.0122; 700.5446  701.0476  680.5617  694.0632  665.1878; 724.8864  719.3337  720.1775  718.7011  682.1620;...
  730.1003  732.2224  736.7990  710.9673  711.6942; 747.0867  755.5361  728.9984  746.5740  736.7782; 778.8731  756.2156  760.1603  758.3994  749.6740;... 
  805.9923  800.7315  786.7603  788.7258  759.4764];


figure;
subplot(1,2,1);
hold on
meanies1 = mean(data1);
meanies2 = mean(data2);
xs1 = 1:length(meanies1);
xs2 = [1:length(meanies1)]+.03;
es1 = std(data1)/sqrt(size(data1,1));
es2 = std(data2)/sqrt(size(data2,1));
errorbar(xs1, meanies1, es1, 'Marker', 'o');
errorbar(xs2, meanies2, es2, 'Marker', '^', 'LineStyle', '--', 'Color', 'r');
tickSpace = 590:20:690;
xLimit = [0.5,5.5];
xTick = 1:5;
set(gca, 'YLim', [590 690]);
set(gca, 'YTick', tickSpace);
set(gca, 'XLim', xLimit);
set(gca, 'XTick', xTick);
xlabel('Second Factor');
ylabel('y');
legend('Level 1', 'Level 2');
title('Conventional SEM');
hold off;

subplot(1,2,2);
hold on;
es1 = withinSubjErr(data1);
es2 = withinSubjErr(data2);
errorbar(xs1, meanies1, es1, 'Marker', 'o');
errorbar(xs2, meanies2, es2, 'Marker', '^', 'LineStyle', '--', 'Color', 'r');
set(gca, 'YLim', [590 690]);
set(gca, 'YTick', tickSpace);
set(gca, 'XLim', xLimit);
set(gca, 'XTick', xTick);
xlabel('Second Factor');
ylabel('y');
legend('Level 1', 'Level 2');
title('Within subject error');
hold off;


