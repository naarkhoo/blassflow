function y = standardize(x)

% make your data have mean 0 and variance 1.

y = (x - repmat(mean(x,1),size(x,1),1))./repmat(std(x),size(x,1),1);
