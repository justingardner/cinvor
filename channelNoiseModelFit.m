% channelNoiseModelFit.m
%
%      usage: channel = channelNoiseModelFit(instanceMatrix,channel)
%         by: justin gardner
%       date: 09/08/16
%    purpose: fits noise model of van Bergen, Ma, Pratte and Jehee, 2015 18:1728-30
%
function channel = channelNoiseModelFit(instanceMatrix,channel,varargin)

% check arguments
if nargin < 2
  help fitChannelNoiseModel
  return
end

% parse args
getArgs(varargin,{'noiseModelGridSearchOnly=0','noiseModelFitTolerence=1','noiseModelGridSteps=10'});

% get residual from channel model
residualMatrix = instanceMatrix - channel.channelResponse*channel.channelWeights;

% set search options
global nIter;
global lastLogLike;lastLogLike = inf;
searchOptions = optimset('Display','none','MaxIter',1000000,'MaxFunEvals',1000000,'TolFun',noiseModelFitTolerence,'TolX',noiseModelFitTolerence);

% first grid search for best starting parameters

% settings that control over what range to grid search
rhoMin = 0;rhoMax = 0.4;
sigmaMin = 0.2;sigmaMax = 0.6;
tauMin = 0.6;tauMax = 1.4;
numSteps = noiseModelGridSteps;

% make the parameter arrays
rho = rhoMin:(rhoMax-rhoMin)/(numSteps-1):rhoMax;
sigma = sigmaMin:(sigmaMax-sigmaMin)/(numSteps-1):sigmaMax;
tau = tauMin:(tauMax-tauMin)/(numSteps-1):tauMax;
[rho sigma tau] = meshgrid(rho,sigma,tau);
nParamValues = length(rho(:));

% do the grid search
nIter = inf; % this will set no display of fit
disp(sprintf('(channelNoiseModelFit) Grid searching over rho [%0.2f %0.2f] sigma [%0.2f %0.2f] tau [%0.2f %0.2f] with %i steps for best initial parameters',rhoMin,rhoMax,sigmaMin,sigmaMax,tauMin,tauMax,numSteps));
for iParamValue = 1:nParamValues
  % get the param values as an array
  params = getNoiseModelInitParams(channel,rho(iParamValue),sigma(iParamValue),tau(iParamValue));
  % compute model likelihood
  initLike(iParamValue) = noiseModelLogLikelihood(params,residualMatrix,channel);
end

% choose best initial params
[~,bestIndex] = min(initLike);
initParams = getNoiseModelInitParams(channel,rho(bestIndex),sigma(bestIndex),tau(bestIndex));
disp(sprintf('(channelNoiseModelFit) Start params: rho=%0.2f sigma=%0.2f tau=%0.2f: %f (fit tolerence: %0.4f)',rho(bestIndex),sigma(bestIndex),tau(bestIndex),initLike(bestIndex),noiseModelFitTolerence));
nIter = 0;

% call fminsearch to find best parameters
if ~noiseModelGridSearchOnly
  [noiseModel.params bestLike] = fminsearch(@noiseModelLogLikelihood,initParams,searchOptions,residualMatrix,channel);
else
  noiseModel.params = initParams;
  bestLike = min(initLike);
end

% make sense of arguments
[noiseModel.covar noiseModel.rho noiseModel.sigma noiseModel.tau] = channelNoiseModelCovar(noiseModel.params,channel);
dispHeader;
disp(sprintf('(channelNoiseModelFit) Best params: rho=%0.2f sigma=%0.2f tau=%0.2f: %f',noiseModel.rho,noiseModel.sigma,mean(noiseModel.tau),bestLike));
dispHeader;

% return channel
channel.noiseModel = noiseModel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    noiseModelLogLikelihood    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglike = noiseModelLogLikelihood(params,residualMatrix,channel)

% convert params from array into stucture for ease
[covar rho sigma tau] = channelNoiseModelCovar(params,channel);

% compute the log likelihood of the residual data (invert sign for fminsearch)
loglike = -sum(log(mvnpdf(residualMatrix,zeros(1,channel.info.numVoxels),covar)));

global nIter;
nIter = nIter + 1;
if mod(nIter,1000) == 1
  global lastLogLike;
  disp(sprintf('%06i: rho=%0.2f sigma=%0.2f tau=%0.2f: %f (delta: %f)',nIter,rho,sigma,mean(tau),loglike,loglike-lastLogLike));
  lastLogLike = loglike;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getNoiseModelInitParams    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = getNoiseModelInitParams(channel,rho,sigma,tau)

% parameters are tau, followed by sigma and n rhos (n is number of voxels)
params = [rho sigma repmat(tau,1,channel.info.numVoxels)];

