% getCrossValSets.m
%
%      usage: crossVal = getCrossValSets(instances,'nFold=5')
%         by: justin gardner
%       date: 09/07/16
%    purpose: gets train and test set of instances
%
function crossVal = getCrossValSets(instances,varargin)

% check arguments
if ~any(nargin == [1])
  help getCrossvalSets
  return
end

% get arguments
getArgs(varargin,{'nFold=5'});

% do n-fold cross-validation. First need to figure out how many instances
% there are in each class
for iClass = 1:length(instances)
  % compute number of instances per class
  crossVal.nInstances(iClass) = size(instances{iClass},1);
end

% save nFold in return strucutre
crossVal.nFold = nFold;

% display number of instances
if unique(crossVal.nInstances)
  disp(sprintf('(getCrossValSets) Number of instances per class is: %i',unique(crossVal.nInstances)));
  balancedInstances = true;
  % only need to compute the train/test sets for one class and then will later
  % replicate.
  nClass = 1;
else
  disp(sprintf('(getCrossValSets) !!! Number of instances is not the same per class. Average is: %0.2f +- %0.2f std',mean(crossVal.nInstances),std(crossVal.nInstances)));
  balancedInstances = false;
  % need to compute train/test separately for each instance
  nClass = length(crossVal.nInstances);
end

% compute train and test
for iClass = 1:nClass
  % compute what instances are going to be used as train and test on each fold
  nTrialsPerFold = crossVal.nInstances(iClass)/nFold;
  % sanity check
  if nTrialsPerFold < 1
    disp(sprintf('(getCrossvalSets) !!! Number of trials per fold less than one of nFold=%i (total number of instances: %i). Set nFold lower',nFold,crossVal.nInstances(iClass)));
    keyboard
  end
  % nTrialsPerFold will not in general be an integer number, so we 
  % have to make them an integer numbers, so we first will
  % randomy add one trial to each fold so that on average we get the
  % non-integer part - i.e. if we need 5.4 trials, then we need
  % to add 1 to each fold 40% of the time.
  nTrialsEachFold = floor(nTrialsPerFold) + (rand(1,nFold) < repmat(nTrialsPerFold-floor(nTrialsPerFold),1,nFold));
  % now fix nTrialsEachFold so that it sums to exactly the number of instances we have
  while sum(nTrialsEachFold) > crossVal.nInstances(iClass)
    % find the fold with the largest number of instances
    % and subtract one from there
    [~,maxIndex] = max(nTrialsEachFold);
    nTrialsEachFold(maxIndex) = nTrialsEachFold(maxIndex)-1;
  end
  % check for smaller than the desired number
  while sum(nTrialsEachFold) < crossVal.nInstances(iClass)
    % find the fold with the smallest number of instances
    % and add one from there
    [~,minIndex] = min(nTrialsEachFold);
    nTrialsEachFold(minIndex) = nTrialsEachFold(minIndex)+1;
  end
  % get a random permutation for which we will take some
  randpermInstances = randperm(crossVal.nInstances(iClass)); 
  % start trial for the train set is 1
  startNum = 1;
  % and the end num is taken from the cumulative sum from nTrialsEachFold
  endNums = cumsum(nTrialsEachFold);
  for iFold = 1:nFold
    % get the end number for the test set for this fold
    endNum = endNums(iFold);
    % and the number for all trials in the test set
    testNums = startNum:endNum;
    % then the train numbers are everything besides these
    trainNums = setdiff(1:crossVal.nInstances(iClass),testNums);
    % now lookup the trial numbers that we will use
    crossVal.trainInstances{iClass,iFold} = randpermInstances(trainNums);
    crossVal.testInstances{iClass,iFold} = randpermInstances(testNums);
    % update start num
    startNum = endNum+1;
  end
  % save TrialsEachFold in crossVal strucutre
  crossVal.nTrialsEachFold{iClass} = nTrialsEachFold;
end

% if we have balanced instances, then replicate
% the train and test values across the rest of the instances
if balancedInstances
  for iClass = 2:length(crossVal.nInstances)
    crossVal.nTrialsEachFold{iClass} = crossVal.nTrialsEachFold{1};
    for iFold = 1:nFold
      crossVal.trainInstances{iClass,iFold} = crossVal.trainInstances{1,iFold};
      crossVal.testInstances{iClass,iFold} = crossVal.testInstances{1,iFold};
    end
  end
end

