function [X_norm, mu, sigma] = featureNormalize(X)
%FEATURENORMALIZE Normalizes the features in X 
%   FEATURENORMALIZE(X) returns a normalized version of X where
%   the mean value of each feature is 0 and the standard deviation
%   is 1. This is often a good preprocessing step to do when
%   working with learning algorithms.

%old_mean is a vector that contains the old mean value of feature1 (in col.1)and
%feature 2 (in col.2)

mu = mean(X); % 1x2

%Create a matrix and compute the standard deviation of each column-feature.
sigma = std(X); % 1x2

Xcentered = X-mu; % 100x2 (implicit expansion by Matlab...) 

X_norm = Xcentered./sigma;


% ============================================================

end
