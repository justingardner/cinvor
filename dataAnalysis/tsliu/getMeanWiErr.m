function [mu, err]=getMeanWiErr(dat, withinDim)
%this function return the mean and std err of a input data matrix, just a
%little utility to save some typing. Note the first dimension of the input
%is assumed to be the number of observations (e.g., subjects, trials)

n=size(dat);
disp(['Assume ', num2str(n(1)), ' observations in the data.']);

nanInData=any(isnan(dat(:)));
if nanInData
  disp('Detected nan in data, using nanmean and nanstd instead');
  error('not implemented for nans in data');
end


mu=squeeze(mean(dat));
if withinDim==2
  for i=1:n(3)
    err(:,i)=withinSubjErr(squeeze(dat(:,:,i)));
  end
elseif withinDim==3
  for i=1:n(2)
    err(i,:)=withinSubjErr(squeeze(dat(:,i,:)));
  end
  disp('THIS NEEDS TO BE TESTED');
end
