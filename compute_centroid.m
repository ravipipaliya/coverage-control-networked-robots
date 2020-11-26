% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
function Cv = compute_centroid(px,py)

global n Fi ai;

% Compute local Voronoi region Vi
% Truncate Vi to Vi_plus
bs_ext = [0 1 1 0; 0 0 1 1]'; % boundary 
[v,c,p] = VoronoiLimit(px ,py,'bs_ext',bs_ext,'figure', 'off'); % Voronoi edges
% Note the function orders region in CW order and changes robot order

Cv = zeros(2,n);
for i = 1:n
    vind = c{i};
    vx = v(vind,1)';
    vy = v(vind,2)';
    pos_x = find(px==p(i,1));
    pos_y = find(py==p(i,2));

% polyin = polyshape({x},{y});
% [a,b] = centroid(polyin);
% plot(polyin)

% g = @(m,n) exp(-(m+n).^2);
% g=@(x,y)exp(x+y);
% g = @(x,y) eval(subs(q(1)+q(1),q,[x;y]))
% g = @(x,y) sqrt(x.^2+y.^2)
% g = @(x,y) 1/(sigma*sqrt(2*pi))*exp(-((x).^2 + (y).^2) /(2*sigma^2));
% intpoly(g,[0,1,2],[0,1,0])
% MonteCarlo_double(
% intpoly(g,x,y)

% k = kappa*Ai(:,pos_x)
%     f = @(x,y) eval(subs(k,q,[x;y]));
%     intpoly(f,[0,1,2],[0,1,0])

    if pos_x == pos_y
        % Reorder centroids based on robot order
        [Cv(:,pos_x)] = dintegrate(vx,vy,p(i,:),pos_x);
    else
        disp('Mismatch in position found')
    end
end
% scatter(Cv(1,:),Cv(2,:),'r*')