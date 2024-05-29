function [ eigenval, eigenvec, order] = myPCA(X)
%PCA Run principal component analysis on the dataset X
%   [ eigenval, eigenvec, order] = mypca(X) computes eigenvectors of the autocorrelation matrix of X
%   Returns the eigenvectors, the eigenvalues (on diagonal) and the order 
%

% Useful values
[m, n] = size(X); %lines x columns 

mu = mean(X); % 1x2

% Make sure each feature from the data is zero mean
X_centered = X-mu; 
mu_centered = mean(X_centered); %should be zero. 

X_trans = X_centered.'; 

covariance_matrix = (X_trans * X_centered)/m ; 

% ====================== YOUR CODE HERE ======================
%

[V,D]=eig(covariance_matrix);

eigenval= diag(D); % Eigenvalues lie in the diagonal of the covariance matrix

eigenvec = V; %Corresponding eigenvectors

[eigenval, order] = sort(eigenval, 'descend');  %Sorting eigenvalues and returning the order
eigenvec = eigenvec(:, order);


% =========================================================================

end
