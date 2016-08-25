function [mu, err]=getMeanErr(dat)
%this function return the mean and std err of a input data matrix, just a
%little utility to save some typing. Note the first dimension of the input
%is assumed to be the number of observations (e.g., subjects, trials)

n=size(dat,1);
disp(['Assume ', num2str(n), ' observations in the data.']);

nanInData=any(isnan(dat(:)));
if nanInData
  disp('Detected nan in data, using nanmean and nanstd instead');
end

if nanInData
  mu=squeeze(nanmean(dat));
  err=squeeze(nanstd(dat)/sqrt(n));
else
  mu=squeeze(mean(dat));
  err=squeeze(std(dat)/sqrt(n));
end