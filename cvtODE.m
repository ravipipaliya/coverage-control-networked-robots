% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% Function that integrates the ODE model governing the robots' dynamics
% ODE for Estimated centroidal voronoi tessalation
function dz = cvtODE(t,z)
global n K Fi ai li Li amin Tau g w_t psi a tprev;
% ---- The Adaptive Coverage Control Algorithm ---- %
px = z(1:n,1);
py = z(n+1:2*n,1);
w_t = exp(-t);
% w_t = (t<2);
% w_t = ones(1,1);

% Voronoi regions & Centroid calcuation
Cv= compute_centroid(px,py); % Updates Fi,Li,li

% T = delaunayTriangulation(px,py);
% ed = edges(T);
% L = zeros(n,n);
% for ind = 1:length(ed)
%     L(ed(ind,1),ed(ind,2))= -1;
% end
% for i = 1:n
%     L(i,i) = -sum(L(i,:));
% end

for i = 1:n
    s = zeros(9,1);
%     for j = 1:n
%         s = L(i,j)*(ai(:,i) - ai(:,j));
%     end
    
    % Update ai
%     size(gamma.*(Li(:,:,i)*ai(:,i) - li(:,i)))
    ai_pre = -(Fi(:,:,i)*ai(:,i)) - g*(Li(:,:,i)*ai(:,i) - li(:,i)) - psi*s;
%     ai_pre'
    Iproji = ones(9,1);
    Iproji(ai(:,i) > amin) = 0;
    Iproji(ai(:,i) == amin & ai_pre >= 0) = 0;
%     size(Iproji)
    ai(:,i) = ai(:,i) + Tau*(ai_pre - diag(Iproji)*ai_pre);
    
end
disp(ai(1:2,:))
% Apply control input: ui = -K(Cvi - pi)
dx = K(1,1)*(Cv(1,:)' - px);
dy = K(2,2)*(Cv(2,:)' - py);

% dx(abs(dx) > K(1,1)) = 0;
% dy(abs(dy) > K(2,2)) = 0;
% disp(norm(a-ai)) % paramter convergence
% disp(ai)
% disp(datestr(now,'HH:MM:SS.FFF'))
dz = [dx; dy;]; 