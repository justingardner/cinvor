% vonMisesFWHM.m
%
%      usage: y = vonMisesFWHM(kappa)
%         by: dylan cable
%       date: 6/2016
%    purpose: 
%
%
%calculated full width half maximum of von mises distribution
function y = vonMisesFWHM(kappa)

param=[pi,kappa,0,1];

xp= 1:180;
xpOri= 1:180;
yp1=vonMises(param,d2r(xp));
yp1_base0=yp1-min(yp1);
rgHalf1=find(yp1_base0>max(yp1_base0)/2);
y=xpOri(rgHalf1(end))-xpOri(rgHalf1(1));