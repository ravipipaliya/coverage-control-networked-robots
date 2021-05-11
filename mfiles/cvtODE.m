% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% Function that integrates the ODE model governing the robots' dynamics
% ODE for Estimated centroidal voronoi tessalation
function dz = cvtODE(t,z)

global n h K Fi amin Tau g psi a kappa;
global est_pos_err tru_pos_err;

% ------------------------------------------------- %
% ---- The Adaptive Coverage Control Algorithm ---- %
% ------------------------------------------------- %
% Initalize
[px,py,ai,li,Li] = reshape_state(z);
dai = zeros(size(ai));
dli = zeros(size(li));
dLi = zeros(size(Li));

% Voronoi regions & Centroid calcuation
[Cv,Cv_true,L] = compute_centroid(px,py,ai); % Computes Fi as well

% Control law: ui = -K(Cvi - pi)
dx = K(1,1)*(Cv(1,:)' - px);
dy = K(2,2)*(Cv(2,:)' - py);
est_pos_err(round(t/h+1)) = mean(vecnorm([(Cv(1,:)' - px),(Cv(2,:)' - py)]'));
tru_pos_err(round(t/h+1)) = mean(vecnorm([(Cv_true(1,:)' - px),(Cv_true(2,:)' - py)]'));

% Adaption laws for paramter estimate
s = ai*L';
for i = 1:n
    dai_pre = -(Fi(:,:,i)*ai(:,i)) - g*(Li(:,:,i)*ai(:,i) - li(:,i)) - psi*s(:,i);
    Iproji = zeros(9,1);
    Iproji(ai(:,i) + dai_pre*h <= amin) = 1;
%     Iproji(ai(:,i) > amin) = 0;
%     Iproji(ai(:,i) == amin & dai_pre >= 0) = 0;
    dai(:,i) = Tau*(dai_pre - diag(Iproji)*dai_pre);    
end

% Update li and Li
% w_t = exp(-t); 
for i = 1:n
    w_t = norm([dx(i) dy(i)])/norm(K); % Data weighting function
    ki = kappa(px(i),py(i));
    phi_t = ki*a;
    dLi(:,:,i) = w_t*(ki'*ki);
    dli(:,i) = w_t*(phi_t)*ki';
end

% State update
dz = [dx; dy; dai(:); dli(:); dLi(:)]; 

%% Debugging ---------------------------- %
disp(t)
disp(mean(vecnorm(a-ai))) % paramter error
% disp(norm([dx,dy])) % distance to centroid
end