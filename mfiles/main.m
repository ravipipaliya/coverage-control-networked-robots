% MAE 598 Multi robot systems
% Project - Ravi Pipaliya
% Distributed Adaptive coverage control
% 11/22/2020
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
addpath('mfiles\','ODE_Solvers\')
clear all;
rng('default'); % Setting seed for reproducibility
global n K Tau g psi h;
global amin Fi kappa a sens_info_flag;
global est_pos_err tru_pos_err; 

%% Initialize ----%
n = 20; % # of robots
K = 3*eye(2); % Control gain matrix
Tau = eye(9); 
g = 130; % Learning rate <100
psi = 0; % Consensus weight <20
h = 0.01; % ode step size
sens_info_flag = 0; % Set 1 to run the case when sensory info is known

%% Initial Positions ----%
x0 = rand(n,1);  % Initial x
y0 = rand(n,1);  % Initial y

%% Model paramters ----%
amin = 0.1; % minimum weight
a = [100 amin*ones(1,7) 100]'; % True weights
ai = amin * ones(9,n); % parameter estimate with each robot
li = zeros(9,n); 
Li = zeros(9,9,n);
Fi = zeros(9,9,n);

%% Sensory basis function ----%
syms qx qy;
sigma = 0.18; % Gaussian sd
mu = [1 1 1 3 3 3 5 5 5; 1 3 5 1 3 5 1 3 5]/6; % Gaussian means
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

%% Simulation ----%
tspan = 0:h:30;
est_pos_err = zeros(length(tspan),1);
tru_pos_err = zeros(length(tspan),1);

z0 = [x0; y0; ai(:); li(:); Li(:)];  % Initial state
z = ode1(@cvtODE,tspan,z0); % Fixed time step ode

%% Save/Load workspace ----%
% Executed takes long so save workspace to compare iterations 

% save('output/base_z.mat') % Sensory distribution is known 
% save('output/z_130_0_err.mat') % Basic non concensus g = 130; psi = 0
% save('output/z_65_15_err.mat') % Concensus g = 65; psi = 15

% load('output/base_z.mat') % Sensory distribution is known 
% load('output/z_130_0_err.mat') % Basic non concensus g = 130; psi = 0
% load('output/z_65_15_err.mat') % Concensus g = 65; psi = 15

%% Decompose state ----%
par_err = zeros(size(z,1),1);
pxi = zeros(n,size(z,1));
pyi = zeros(n,size(z,1));
for i = 1:size(z,1)
    [pxi(:,i),pyi(:,i),ain] = reshape_state(z(i,:)');
    par_err(i) = mean(vecnorm(a-ain));
end

%% Plots ----%
% a. Initial configuration
figure
subplot(2,2,1)
voronoi(x0,y0,'b.')
hold on
plot(mu(1,1),mu(2,1),'r*')
plot(mu(1,9),mu(2,9),'r*')
title('Initial configuration')

% b. Final configuration
subplot(2,2,2)
voronoi(pxi(:,end),pyi(:,end),'b.')
hold on
plot(mu(1,1),mu(2,1),'r*')
plot(mu(1,9),mu(2,9),'r*')
title('Final configuration')

% c. Parameter convergence analysis
subplot(2,2,3)
plot(tspan,par_err,'-r')
title('Parameter convergence')
ylabel('$$Mean ||\tilde{a}_i(t)||$$','interpreter','latex')
xlabel('time')
ylim([0,150])

% d. Robot trajectories
subplot(2,2,4)
hold on;
for i=1:n
    plot(pxi(i,end)',pyi(i,end)','ko')
    plot(pxi(i,:)',pyi(i,:)','--')
end
title('Robot trajectories')

%% Gifs ----%
% Defining the range of axes
% ax = axes('XLim',[0 1],'YLim',[0 1]);
cf = gcf;
filename = 'gifs/consensus_test.gif';
for i = 1:size(pxi,2)
    hold on
    voronoi(pxi(:,i),pyi(:,i),'b.')
    plot(mu(1,1),mu(2,1),'r*')
    plot(mu(1,9),mu(2,9),'r*')
    title(strcat('t = ',num2str(tspan(i))))
    xlabel('Consensus controller')
    xlim([0,1])
    ylim([0,1])
    hold off
    drawnow
    % Capture the plot as an image 
    frame = getframe(cf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1
        imwrite(imind,cm,filename, 'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.0001,'WriteMode','append'); 
    end
    clf;
end

%% Comparing Position error (Basic non concesus vs concensus controller)
load('output/z_130_0_err.mat') % non concensus g = 130; psi = 0
est_pos_err_nc = est_pos_err;
tru_pos_err_nc = tru_pos_err;
load('output/z_65_15_err.mat') % concensus g = 65; psi = 15
tx = tspan(:,1:end-1);
ty = 1:size(tx,2);

% Estimated position error
figure
subplot(1,2,1)
plot(tx,est_pos_err_nc(ty), 'k')
hold on;
plot(tx,est_pos_err(ty),'r')
title('Mean Estimated Position Error')
legend('Basic', 'Consensus')
ylabel('$$||\hat{C}_{V_i} - p_i||$$','interpreter','latex')
xlabel('time')

% True position error
subplot(1,2,2)
plot(tx,tru_pos_err_nc(ty),'k')
hold on
plot(tx,tru_pos_err(ty), 'r')
title('Mean True Position Error')
legend('Basic', 'Consensus')
ylabel('$$||C_{V_i} - p_i||$$','interpreter','latex')
xlabel('time')