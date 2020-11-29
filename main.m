% MAE 598 Multi robot systems
% Project - Ravi Pipaliya
% Distributed Adaptive coverage control
% 11/22/2020
% clear all;
global n K Tau g psi h;
global amin Fi kappa a;

%---- Initialize ----%
n = 20; % # of robots
K = 3*eye(2); % Control gain matrix
Tau = eye(9); 
g = 65; % Learning rate <100
psi = 15; % Consensus weight <20
h=0.01; % ode step size

%---- Initial Positions ----%
rng('default'); % Setting seed for reproducibility
x0 = 1*rand(n,1);  % Initial x
y0 = 1*rand(n,1);  % Initial y

%---- Model paramters ----%
amin = 0.1; % minimum weight
a = [100 amin*ones(1,7) 100]'; % True weights
ai = amin * ones(9,n); % 9X1 for 
li = zeros(9,n);
Li = zeros(9,9,n);
Fi = zeros(9,9,n);

%---- Basis function ----%
% q = sym('q',[2,1]);
syms qx qy;
sigma = 0.18; % Gaussian sd
mu = [1 1 1 3 3 3 5 5 5; 
      1 3 5 1 3 5 1 3 5]/6;
% gd = @(m)1/(sigma*sqrt(2*pi))*exp(-((qx-m(1)).^2 + (qy-m(2)).^2) /(2*sigma^2));
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

%---- Simulation ----%
tspan = 0:h:10;
z0 = [x0; y0; ai(:); li(:); Li(:)];  % Initial state
% [t,z] = ode45(@cvtODE,tspan,z0);
z = ode1(@cvtODE,tspan,z0);

%---- Plots ----%
% Initial configuration
figure
voronoi(x0,y0,'b.')
hold on
plot(mu(1,1),mu(2,1),'r*')
plot(mu(1,9),mu(2,9),'r*')
xlim([0 1])
ylim([0 1])

% Final configuration
xn = z(end,1:n);
yn = z(end,n+1:2*n);
[pxn,pyn,ain,lin,Lin] = reshape_state(z(end,:)');

figure
voronoi(pxn,pyn,'b.')
hold on
plot(mu(1,1),mu(2,1),'r*')
plot(mu(1,9),mu(2,9),'r*')
xlim([0 1])
ylim([0 1])

% Parameter analysis
err = zeros(size(z,1),1);
pxi = zeros(n,size(z,1));
pyi = zeros(n,size(z,1));
for i = 1:size(z,1)
    [pxi(:,i),pyi(:,i),ain] = reshape_state(z(i,:)');
    err(i) = mean(vecnorm(a-ain));
end
figure;
plot(tspan,err) % Paramter convergence

% Robot trajectories
figure;
voronoi(pxn,pyn,'b.');
hold on;
for i=1:n
    plot(pxi(i,1)',pyi(i,1)','ko')
    plot(pxi(i,:)',pyi(i,:)','--')
end