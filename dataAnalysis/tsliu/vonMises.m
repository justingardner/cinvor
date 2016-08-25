% vonMises.m
%
%      usage: y = vonMises(params, x)
%         by: taosheng liu
%       date: 11/2015
%    purpose: 
%
%
%a 4-parameter circular normal (von Mises) function
%params(1):mean; params(2):kappa; params(3):baseline; params(4):amplitude
function y = vonMises(params, x)

% check arguments
if ~any(nargin == [2])
  help vonMises
  return
end

% decode parameters
mu=params(1);
kappa=params(2);
base=params(3);
amp=params(4);

% use libary function to calculate von mises
vm=circ_vmpdf(x,mu,kappa);

% normalize
vm_norm=(vm-min(vm))/(max(vm)-min(vm));

% and on baseline and amplitude
y=base+amp*vm_norm;

%circ_vmpdf always return column vector, turn it into row
y=y'; 
