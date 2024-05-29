%TEL311
close all;
clear;
clc;

%Exercise 1 - Bayes
% Parameters of the two Gaussian functions
mu1 = [2; 3];
mu1
sigma1 = [2 0.5; 0.5 1];
sigma1
mu2 = [4; 4];
mu2
sigma2 = [1.5 -0.3; -0.3 0.8];
sigma2

% Define the grid of points to evaluate the Gaussian functions and the decision boundary
x = linspace(-5, 10, 100);
y = linspace(-5, 10, 100);
[X, Y] = meshgrid(x, y);

% Evaluate the Gaussian functions at each point on the grid
z1 = mvnpdf([X(:) Y(:)], mu1', sigma1);
Z1 = reshape(z1, size(X));

z2 = mvnpdf([X(:) Y(:)], mu2', sigma2);
Z2 = reshape(z2, size(X));

%plot only the distributions : 
figure() 
% Plot the 3D Gaussian functions
surf(X, Y, Z1, 'FaceColor', 'red');
hold on;
surf(X, Y, Z2, 'FaceColor', 'blue');
grid on;
xlabel('x');
ylabel('y');
zlabel('Probability density');
title('Two 3D Gaussian Functions');
hold off;

% Compute the decision boundary for different values of P(?1)
P_omega1 = [0.1, 0.25, 0.5, 0.75, 0.9];
colors = {'g', 'm', 'c', 'y', 'k'};
figure()
hold on;
for i = 1:length(P_omega1)
    log_term = log(P_omega1(i)/(1-P_omega1(i)));
    Z = log(Z1) - log(Z2) + log_term;
    contour(X, Y, Z, [0 0], 'LineWidth', 2, 'LineColor', colors{i});
end

% Plot the 3D Gaussian functions
surf(X, Y, Z1, 'FaceColor', 'red');
surf(X, Y, Z2, 'FaceColor', 'blue');
grid on;
xlabel('x');
ylabel('y');
zlabel('Probability density');
title('Two 3D Gaussian Functions (Case 1)');
legend('P1=0.1', 'P1=0.25', 'P1=0.5', 'P1=0.75', 'P1=0.9', 'Gaussian 1', 'Gaussian 2');
hold off;

%==============================================================
fprintf('\nCase 2\n\n');

mu1 = [2; 3];
mu1
sigma1 = [1.2 0.4; 0.4 1.2];
sigma1
mu2 = [4; 4];
mu2
sigma2 = [1.2 0.4; 0.4 1.2];
sigma2


% Define the grid of points to evaluate the Gaussian functions and the decision boundary
x = linspace(-5, 10, 100);
y = linspace(-5, 10, 100);
[X, Y] = meshgrid(x, y);

% Evaluate the Gaussian functions at each point on the grid
z1 = mvnpdf([X(:) Y(:)], mu1', sigma1);
Z1 = reshape(z1, size(X));

z2 = mvnpdf([X(:) Y(:)], mu2', sigma2);
Z2 = reshape(z2, size(X));

%plot only the distributions : 
figure() 
% Plot the 3D Gaussian functions
surf(X, Y, Z1, 'FaceColor', 'red');
hold on;
surf(X, Y, Z2, 'FaceColor', 'blue');
grid on;
xlabel('x');
ylabel('y');
zlabel('Probability density');
title('Two 3D Gaussian Functions');
hold off;

% Compute the decision boundary for different values of P(?1)
P_omega1 = [0.1, 0.25, 0.5, 0.75, 0.9];
colors = {'g', 'm', 'c', 'y', 'k'};
figure()
hold on;
for i = 1:length(P_omega1)
    log_term = log(P_omega1(i)/(1-P_omega1(i)));
    Z = log(Z1) - log(Z2) + log_term;
    contour(X, Y, Z, [0 0], 'LineWidth', 2, 'LineColor', colors{i});
end

% Plot the 3D Gaussian functions
surf(X, Y, Z1, 'FaceColor', 'red');
surf(X, Y, Z2, 'FaceColor', 'blue');
grid on;
xlabel('x');
ylabel('y');
zlabel('Probability density');
title('Two 3D Gaussian Functions (Case 2)');
legend('P1=0.1', 'P1=0.25', 'P1=0.5', 'P1=0.75', 'P1=0.9', 'Gaussian 1', 'Gaussian 2');
hold off;