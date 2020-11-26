% MAE 598 Multi robot systems
% Project - Ravi Pipaliya
% Distributed Adaptive coverage control
% 11/22/2020
clear all;
global n K amin Fi ai li Li kappa a Tau g psi;

%---- Initialize ----%
n = 4; % # of robots
K = 3*eye(2); % Control gain matrix
Tau = 0.01*eye(9); 
g = 300; % Learning rate
psi = 0; % Consensus weight

%---- States ----%
rng('default'); % Setting seed for reproducibility
x0 = 1*rand(n,1);  % Initial x
y0 = 1*rand(n,1);  % Initial y
z0 = [x0; y0;];  % Initial state

%---- Model paramters ----%
amin = 0.1; % minimum weight
a = [100 amin*ones(1,7) 100]'; % True weights
ai = amin * ones(9,n); % Inital parameters
li = zeros(9,n);
Li = zeros(9,9,n);
Fi = zeros(9,9,n);

%---- Basis function ----%
% q = sym('q',[2,1]);
syms qx qy;
sigma = 0.18; % Gaussian sd
% gd = @(m)1/(sigma*sqrt(2*pi))*exp(-((qx-m(1)).^2 + (qy-m(2)).^2) /(2*sigma^2));
mu = [1 1 1 3 3 3 5 5 5; 
      1 3 5 1 3 5 1 3 5]/6;
% kappa = splitapply(gd,mu,1:size(mu,2))';
kappa = @(qx,qy) ...
   [1/(sigma^2*(2*pi))*exp(-((qx-mu(1,1)).^2 + (qy-mu(2,1)).^2) /(2*sigma^2)),...
    1/(sigma^2*(2*pi))*exp(-((qx-mu(1,2)).^2 + (qy-mu(2,2)).^2) /(2*sigma^2)),...
    1/(sigma^2*(2*pi))*exp(-((qx-mu(1,3)).^2 + (qy-mu(2,3)).^2) /(2*sigma^2)),...
    1/(sigma^2*(2*pi))*exp(-((qx-mu(1,4)).^2 + (qy-mu(2,4)).^2) /(2*sigma^2)),...
    1/(sigma^2*(2*pi))*exp(-((qx-mu(1,5)).^2 + (qy-mu(2,5)).^2) /(2*sigma^2)),...
    1/(sigma^2*(2*pi))*exp(-((qx-mu(1,6)).^2 + (qy-mu(2,6)).^2) /(2*sigma^2)),...
    1/(sigma^2*(2*pi))*exp(-((qx-mu(1,7)).^2 + (qy-mu(2,7)).^2) /(2*sigma^2)),...
    1/(sigma^2*(2*pi))*exp(-((qx-mu(1,8)).^2 + (qy-mu(2,8)).^2) /(2*sigma^2)),...
    1/(sigma^2*(2*pi))*exp(-((qx-mu(1,9)).^2 + (qy-mu(2,9)).^2) /(2*sigma^2))];
% gd = @(x,y) kappa(x,y)*a
% fsurf(gd,[0,1])
%---- Simulation ----%
[t,z] = ode45(@cvtODE,[0 1],z0);

%---- Plots ----%
% Initial configuration
figure
voronoi(x0,y0,'b.')
hold on
plot(mu1(1),mu1(2),'r*')
plot(mu2(1),mu2(2),'r*')
xlim([0 1])
ylim([0 1])

% Final configuration
xn = z(end,1:n);
yn = z(end,n+1:2*n);
figure
voronoi(xn,yn,'b.')
hold on
plot(mu1(1),mu1(2),'r*')
plot(mu2(1),mu2(2),'r*')
xlim([0 1])
ylim([0 1])