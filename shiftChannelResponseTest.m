% shift.m
%
%      usage: shift()
%         by: justin gardner
%       date: 04/25/17
%    purpose: 
%
function retval = shift()


% compute filters over the following x
x = d2r(0:0.1:180);

% desired stimulus value
stimVal1 = d2r(90);
[~,stimIndex1] = min(abs(x - stimVal1));

% actual stimulus value
stimVal2 = d2r(78);
[~,stimIndex2] = min(abs(x - stimVal2));

% basis filters for original and desired shift
exponent = 7;
[fun filterPref] = makeFilters(x,exponent,0);
nFilters = length(filterPref);

% compute response to stimulus
funResponse1 = fun(stimIndex1,:);
funResponseScaled1 = fun.*repmat(funResponse1,length(x),1);
funResponse2 = fun(stimIndex2,:);
funResponseScaled2 = fun.*repmat(funResponse2,length(x),1);

% where to recenter the response
recenterPos = stimVal1;

% make the preference difference matrix which is the difference in prefered filter preference
% from the desired position
prefDiffMatrix = mod(repmat(filterPref',1,nFilters) - repmat(filterPref,nFilters,1)-(stimVal1-stimVal2),pi);
prefDiffMatrix(prefDiffMatrix>pi/2) = pi-prefDiffMatrix(prefDiffMatrix>pi/2);

% compute a weight matrix which gives the linear weights to the two
% nearest neighbor filters. (e.g. should be 0.5 and 0.5 if the new
% location is inbetween the two).
filterSpacing = median(diff(filterPref));
weightMatrix = 1 - abs(prefDiffMatrix)/filterSpacing;
weightMatrix(abs(prefDiffMatrix)>=filterSpacing) = 0;

% now get the expected filter response for the original location
% and the shifted location
shiftedFilterResponse = interp1(x,fun,recenterPos);
originalFilterResponse =  interp1(x,fun,stimVal2);
% the ratio is then the expectation for how the new locations response
% should depend on the original response
interpMatrix = repmat(shiftedFilterResponse,nFilters,1) ./ repmat(originalFilterResponse',1,nFilters);

% this interpation blows up whenever the filter vaules go to zero, so take those out
weightMatrix(isinf(interpMatrix)) = 0;
weightMatrix(isnan(interpMatrix)) = 0;
interpMatrix(isnan(interpMatrix)) = 0;
interpMatrix(isinf(interpMatrix)) = 0;

% this may take out some values that should have been weighted, so need to 
% fix that so that the weights all sum to 1 in each column
sumOfWeights = sum(weightMatrix);
sumOfWeights(sumOfWeights==0) = 1;
weightMatrix = weightMatrix * 1./repmat(sumOfWeights,nFilters,1);

% the shift matrix is now just the linear weights of neigbors
% by how much each shoudl be scaled by
shiftMatrix = weightMatrix.*interpMatrix

% display
mlrSmartfig('shifttest','reuse');clf;

% display dimensions
xMax = 180;
yMin = 0;yMax = 1;

% the response we are trying to recreate
subplot(3,1,1);
plot(r2d(x),funResponseScaled1);hold on
yaxis(yMin,yMax);
vline(r2d(stimVal1));
title('Desired interpolation');

% the response we got for the stimulus
subplot(3,1,2);
plot(r2d(x),funResponseScaled2);
hold on
yaxis(yMin,yMax);
vline(r2d(stimVal2));
title('New stimulus');
ylabel('Filter response');

% the interpolated response
subplot(3,1,3);
plot(r2d(x),fun.*repmat(funResponse2*shiftMatrix,length(x),1));
yaxis(yMin,yMax);
title('Interpolated response');
xlabel('Orientation');

keyboard

%%%%%%%%%%%%%%%%%%%%%
%    makeFilters    %
%%%%%%%%%%%%%%%%%%%%%
function [fun filterPref] = makeFilters(x,exponent,phaseShift)

fun = [];
phaseVals = (0:2*pi/((exponent+1)):(2*pi-0.001))+phaseShift;

for iPhase = 1:length(phaseVals)
  % cos function
  fun(:,iPhase) = cos(2*x-phaseVals(iPhase));
  % rectify
  fun(fun(:,iPhase)<0,iPhase) = 0;
  % exponent
  fun(:,iPhase) = fun(:,iPhase).^exponent;
end

filterPref = phaseVals/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeFiltersNonRect    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = makeFiltersNonRect(x,exponent,phaseShift)

fun = [];
phaseVals = (0:pi/(exponent+1):(pi-0.001))+phaseShift;

for iPhase = 1:length(phaseVals)
  % cos function
  fun(:,iPhase) = cos(x+phaseVals(iPhase));
  % exponent
  fun(:,iPhase) = fun(:,iPhase).^exponent;
end
