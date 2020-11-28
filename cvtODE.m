% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% Function that integrates the ODE model governing the robots' dynamics
% ODE for Estimated centroidal voronoi tessalation
function dz = cvtODE(t,z)
global n h K Fi amin Tau g psi a kappa;

% ------------------------------------------------- %
% ---- The Adaptive Coverage Control Algorithm ---- %
% ------------------------------------------------- %
[px,py,ai,li,Li] = reshape_state(z);
% ai(ai<amin) = amin;
dai = zeros(size(ai));
dli = zeros(size(li));
dLi = zeros(size(Li));

% ---- Voronoi regions & Centroid calcuation ---- %
Cv = compute_centroid(px,py,ai); % Computes Fi as well
% voronoi(px,py,'b.');
% hold on;
% scatter(Cv(1,:),Cv(2,:));

% ---- Apply control input: ui = -K(Cvi - pi)- %
dx = K(1,1)*(Cv(1,:)' - px);
dy = K(2,2)*(Cv(2,:)' - py);
% dx(abs(dx) > K(1,1)) = 0;
% dy(abs(dy) > K(2,2)) = 0;

% ---- Update ai ------------------------------- %
T = delaunayTriangulation(px,py);
ed = edges(T);
L = zeros(n,n);
for ind = 1:length(ed)
    L(ed(ind,1),ed(ind,2))= -1;
    L(ed(ind,2),ed(ind,1))= -1;
end
L = L + diag(-1*sum(L,2));
s = ai*L';
% disp(ai(1,1))
% disp(norm(s))

for i = 1:n
    dai_pre = -(Fi(:,:,i)*ai(:,i)) - g*(Li(:,:,i)*ai(:,i) - li(:,i)) - psi*s(:,i);
    Iproji = zeros(9,1);
%     Iproji(dai_pre <= 0 & ai(:,i) < amin) = 1;
%     Iproji(dai_pre <= 0) = 1;
    Iproji(ai(:,i) + dai_pre*h <= amin) = 1;
%     Iproji(ai(:,i) > amin) = 0;
%     Iproji(ai(:,i) == amin & dai_pre >= 0) = 0;
    dai(:,i) = Tau*(dai_pre - diag(Iproji)*dai_pre);    
end

% ---- Update li and Li ---------------------- %
% w_t = exp(-t); % Data weighting function
% w_t = (t<4);
% w_t = ones(1,1);

for i = 1:n
    w_t = norm([dx(i) dy(i)]);
    ki = kappa(px(i),py(i));
    phi_t = ki*a;
    dLi(:,:,i) = w_t*(ki'*ki);
    dli(:,i) = w_t*(phi_t)*ki';
end

dz = [dx; dy; dai(:); dli(:); dLi(:)]; 

% ---- Debugging ---------------------------- %
% disp(px(1))
disp(t)
% disp(ai(1:2,:))
disp(mean(vecnorm(a-ai))) % paramter convergence
% disp(ai)
% disp(datestr(now,'HH:MM:SS.FFF'))
