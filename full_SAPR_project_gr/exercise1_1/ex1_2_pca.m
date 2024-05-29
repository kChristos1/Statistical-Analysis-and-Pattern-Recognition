%% Pattern Recognition
%  Exercise | Principle Component Analysis
%
%  Instructions
%  ------------
%
%  This file contains code that helps you get started on the
%  exercise. You will need to complete the following functions:
%
%     myPCA.m
%     projectData.m
%     recoverData.m
%
%  Add code where needed.
%

%% Initialization
clear all; close all; 

%% ================== Part 1: Load Example Dataset  ===================
%  We start this exercise by using a small dataset that is easily to
%  visualize
%
fprintf('Visualizing example dataset for PCA.\n\n');


%% a1
%  The following command loads the breast cancert dataset. You should now have the 
%  variable X in your environment
Data=csvread('breast_cancer_data.csv');
 
NSamples = 300; %Get the first NSamples only for better visualization

X=Data(1:NSamples,1:end-1); % Get all features

%Y stores values (0/1) for bad or good sample.
Y=Data(1:NSamples,end);

% Visualize the example dataset in 2D using sample feature pairs
figure;

% blue is good 
% red is bad 
plot(X(Y==0, 1), X(Y==0, 2), 'bo',X(Y==1, 1), X(Y==1, 2), 'ro' );
axis square;
title('Data Samples in 2D')
xlabel('Feature 1');
ylabel('Feature 2');
hold on

fprintf('Program paused. Press enter to continue.\n');
pause;

%% b1
% Repeat using different feature pairs on the 2D space 
figure ; 
subplot(1,2,1)
plot(X(Y==0, 3), X(Y==0, 4), 'bo',X(Y==1, 3), X(Y==1, 4), 'ro' );
axis square;
title('Data Samples in 2D')
xlabel('Feature 3');
ylabel('Feature 4');
hold on
subplot(1,2,2)
plot(X(Y==0, 27), X(Y==0, 28), 'bo',X(Y==1, 27), X(Y==1, 28), 'ro' );
axis square;
title('Data Samples in 2D')
xlabel('Feature 27');
ylabel('Feature 28');
hold on

%% a2
%  You should now implement PCA, a feature extraction algorithm. 
%  Complete the code in myPCA.m
%  Consider only 2D samples by taking the first two features of the dataset

fprintf('\nRunning PCA on example dataset taking the first 2 features.\n\n');
X=Data(1:NSamples,1:2); % Get 2 first features

% Before running PCA, it is important to first normalize X
% Add the necessary code to perform standardization in featureNormalize
[X_norm, mu, sigma] = featureNormalize(X);

%plotting the normalized samples
figure ; 
% blue is good 
% red is bad 
plot(X_norm(Y==0, 1), X_norm(Y==0, 2), 'bo',X_norm(Y==1, 1), X_norm(Y==1, 2), 'ro' );
axis square;
title('Normalized Data Samples in 2D')
xlabel('Feature 1');
ylabel('Feature 2');



%% b2 
%plotting the normalized samples
figure ; 
%hold on
% blue is good 
% red is bad 
plot(X_norm(Y==0, 1), X_norm(Y==0, 2), 'bo',X_norm(Y==1, 1), X_norm(Y==1, 2), 'ro' );
hold on
axis square;
title('Normalized Data Samples in 2D')
xlabel('Feature 1');
ylabel('Feature 2');
[eigvals, eigvecs, order] = myPCA(X_norm);
% Plot 2 of the principal component vectors (arrows)
scale = 3; % Scale factor for the arrows
quiver(0, 0, eigvecs(1,1)*scale*eigvals(1), eigvecs(2,1)*scale*eigvals(1), 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
quiver(0, 0, eigvecs(1,2)*scale*eigvals(2), eigvecs(2,2)*scale*eigvals(2), 'b', 'LineWidth', 2, 'MaxHeadSize', 1);
title('Normalized Samples and principal components')
xlabel('Normalized feature 1')
ylabel('Normalized feature 2')
grid on;
hold off


%% c, d (use both components) 

% Calculate the Explained Variance
% ADD YOUR CODE
ExplainedVar = eigvals / sum(eigvals);
fprintf(' Explained Variance(1st PC = %f) (2nd PC = %f)\n', ExplainedVar(1), ExplainedVar(2));

% dimention reduction using pca..
% Projection of the data onto the principal components
% ADD YOUR CODE
K=2;
X_PCA = projectData(X_norm,eigvecs,K);

% Plot the projection of the data onto the principal components
figure;%5
hold on
plot(X_PCA(Y==0, 1), X_PCA(Y==0, 2), 'bo','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerSize', 4 )
plot(X_PCA(Y==1, 1), X_PCA(Y==1, 2), 'ro','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerSize', 4 );
title('PCA of Breast Cancer Dataset');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
axis square;
hold off;

%%d
%  Use the K principal components to recover the Data in 2D
X_recover  = recoverData(X_PCA, eigvecs, K);
fprintf('Approximation of the first example: %f %f\n', X_recover(1, 1), X_recover(1, 2));

%  Draw the lines connecting the projected points to the original points
figure;%6
hold on;
plot(X_recover(Y==0, 1), X_recover(Y==0, 2), 'co','MarkerFaceColor', 'c', 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
plot(X_recover(Y==1, 1), X_recover(Y==1, 2), 'yo', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none', 'MarkerSize', 7);

plot(X_norm(Y==0, 1), X_norm(Y==0, 2), 'bo','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerSize', 4 )
plot(X_norm(Y==1, 1), X_norm(Y==1, 2), 'ro','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerSize', 4 );
title('PCA sample projections (2 principal components)');
xlabel('Normalized Feature 1');
ylabel('Normalized Feature 2');

axis square
for i = 1:size(X_norm, 1)
    drawLine(X_norm(i,:), X_recover(i,:), '--k', 'LineWidth', 1);
end
axis([-3 3.5 -3 3.5]); 
axis square
grid on;
hold off

% Print the first principal component
fprintf(' 1st Principal Component = %f %f \n', eigvecs(1,1), eigvecs(2,1));
fprintf(' 2nd Principal Component = %f %f \n', eigvecs(1,2), eigvecs(2,2));
fprintf('Program paused. Press enter to continue.\n');
pause; 



%% use only the first component 
%% =================== Part 3: Dimensionality Reduction ===================
%  Perform dimensionality reduction by projecting the data onto the 
%  first principal component. Then plot the data in this reduced in 
%  dimension space.  This will show you what the data looks like when 
%  using only the corresponding eigenvectors to reconstruct it.
%
%  You should complete the code in projectData.m
%
fprintf('\nDimensionality reduction on example dataset.\n\n');

%  Project the data onto K = 1 dimension
K = 1;
Z = projectData(X_norm, eigvecs, K);
fprintf('Projection of the first example: %f\n', Z(1));

%  Use the K principal components to recover the Data in 2D
X_recover  = recoverData(Z, eigvecs, K);
fprintf('Approximation of the first example: %f %f\n', X_recover(1, 1), X_recover(1, 2));

%  Draw the lines connecting the projected points to the original points
figure;
hold on;
plot(X_recover(Y==0, 1), X_recover(Y==0, 2), 'co','MarkerFaceColor', 'c', 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
plot(X_recover(Y==1, 1), X_recover(Y==1, 2), 'yo', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none', 'MarkerSize', 5);

plot(X_norm(Y==0, 1), X_norm(Y==0, 2), 'bo','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerSize', 5 )
plot(X_norm(Y==1, 1), X_norm(Y==1, 2), 'ro','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerSize', 5 );
title('PCA sample projections');
xlabel('Normalized Feature 1');
ylabel('Normalized Feature 2');

axis square
for i = 1:size(X_norm, 1)
    drawLine(X_norm(i,:), X_recover(i,:), '--k', 'LineWidth', 1);
end
axis([-3 3.5 -3 3.5]); 
axis square
grid on;
hold off

fprintf('Program paused. Press enter to continue.\n');
pause;


%% =============== Part 4: Perform PCA on all features in the dataset =============
%  Apply PCA on the initial dataset and then choose the first 2 Principal Componets
%  to reduce the dimensionality of the dataset to the 2D space
%  The following code will load the dataset into your environment
%

X=Data(1:NSamples,1:end-1); % Get all features from the breast cancer dataset

% Before running PCA, it is important to first normalize X
% Add the necessary code to perform standardization in featureNormalize
[X_norm, mu, sigma] = featureNormalize(X);

%  Run PCA
[eigvals, eigvecs, order] = myPCA(X_norm);

figure %8
hold on
% Plot the samples on the first two principal components
plot(X_norm(Y==0, order(1)), X_norm(Y==0, order(2)), 'bo','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerSize', 5 )
plot(X_norm(Y==1, order(1)), X_norm(Y==1, order(2)), 'ro','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerSize', 5 );
title('PCA sample projections');
xlabel('Normalized Feature 1');
ylabel('Normalized Feature 2');
% Plot the principal component vectors (arrows)
scale = 1; % Scale factor for the arrows
quiver(0, 0, eigvecs(1,1)*scale*eigvals(1), eigvecs(2,1)*scale*eigvals(1), 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
quiver(0, 0, eigvecs(1,2)*scale*eigvals(2), eigvecs(2,2)*scale*eigvals(2), 'b', 'LineWidth', 2, 'MaxHeadSize', 1);
hold off;
axis square;

% Projection of the data onto the principal components
% ADD YOUR CODE
K=2;
X_PCA = projectData(X_norm,eigvecs,K);



% Plot the samples on the first two principal components
figure;%9
scatter(X_PCA(:,1), X_PCA(:,2), 30, Y, 'filled');
title('Breast Cancer Dataset - PCA reduced in 2D');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
hold off;

fprintf(' 1st Principal Component = %f %f \n', eigvecs(1,1), eigvecs(2,1));
fprintf(' 2nd Principal Component = %f %f \n', eigvecs(1,2), eigvecs(2,2));

% Calculate the Explained Variance of the first 2 Principal Components
% ADD YOUR CODE
ExplainedVar = eigvals / sum(eigvals);
fprintf(' Explained Variance(1st PC = %f) (2nd PC = %f)\n', ExplainedVar(1), ExplainedVar(2));
% Print the first principal component
fprintf(' 1st Principal Component = %f %f \n', eigvecs(1,1), eigvecs(2,1));
fprintf('Program paused. Press enter to continue.\n');
pause;


%% =============== Part 5: Loading and Visualizing Face Data =============
%  We start the exercise by first loading and visualizing the dataset.
%  The following code will load the dataset into your environment
%
fprintf('\nLoading face dataset.\n\n');

%  Load Face dataset
load ('faces.mat')


%  Display the first 100 faces in the dataset
figure; 
displayData(X(1:100, :));

fprintf('Program paused. Press enter to continue.\n');
pause;

%% =========== Part 6: PCA on Face Data: Eigenfaces  ===================
%  Run PCA and visualize the eigenvectors which are in this case eigenfaces
%  We display the first 36 eigenfaces.
%
fprintf(['\nRunning PCA on face dataset.\n' ...
         '(this mght take a minute or two ...)\n\n']);

%  Before running PCA, it is important to first normalize X by subtracting 
%  the mean value from each feature
[X_norm, mu, sigma] = featureNormalize(X);

%  Run PCA
[eigvals, eigvecs, order] = myPCA(X_norm);

%  Visualize the top 36 eigenvectors found (eigenfaces)
figure; 
displayData(eigvecs(:, 1:36)');

fprintf('Program paused. Press enter to continue.\n');
pause;


%% ============= Part 7: Dimensionality Reduction for Faces =================
%  Project images to the eigen space using the top k eigenvectors 
%  If you are applying a machine learning algorithm 
fprintf('\nDimension reduction for face dataset.\n\n');

K = 200;
Z = projectData(X_norm, eigvecs, K);

fprintf('The projected data Z has a size of: ')
fprintf('%d ', size(Z));

%fprintf('\n\nProgram paused. Press enter to continue.\n');
%pause;

%% ==== Part 8: Visualization of Faces after PCA Dimensionality Reduction ====
%  Project images to the eigen space using the top K eigen vectors and 
%  visualize only using those K dimensions
%  Compare to the original input, which is also displayed

fprintf('\nVisualizing the projected (reduced dimension) faces.\n\n');

K = 200;

X_rec  = recoverData(Z, eigvecs, K);

% Display normalized data
figure;
subplot(1, 2, 1);
displayData(X_norm(1:100,:));
title('Original faces');
axis square;

% Display reconstructed data from only k eigenfaces
subplot(1, 2, 2);
displayData(X_rec(1:100,:));
title('Recovered faces');
axis square;




