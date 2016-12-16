% channelNoiseModelCovar.m
%
%      usage: [covar rho sigma tau] = channelNoiseModelCovar(params,channel)
%         by: justin gardner
%       date: 09/08/16
%    purpose: Helper function for channelNoiseModelFit which returns the covariance
%             matrix from params
%
function [covar rho sigma tau] = channelNoiseModelCovar(params,channel)

% check arguments
if ~any(nargin == [2])
  help channelNoiseModelCovar
  return
end

% interpert parameters
rho = params(1);
sigma = params(2);
tau = params(3:end)';

% make the covariance matrix
covar = rho*tau*tau' + (1-rho)*diag(tau.*tau) + sigma^2 * channel.channelWeights' * channel.channelWeights;

