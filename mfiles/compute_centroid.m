% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %

function [Cv,Cv_true,L] = compute_centroid(px,py,ai)

global n;

% Compute local Voronoi region Vi
bs_ext = [0 1 1 0; 0 0 1 1]'; % Environment boundary 
[v,c,p] = VoronoiLimit(px ,py,'bs_ext',bs_ext,'figure', 'off'); 
% Note the function orders region in CW order and changes robot order

Cv = zeros(2,n);
Cv_true = zeros(2,n);
order = zeros(n,1);
for i = 1:n
    vind = c{i};
    vx = v(vind,1)';
    vy = v(vind,2)';
    pos_x = find(px==p(i,1));
    pos_y = find(py==p(i,2));
    order(pos_x) = i;
    
    if pos_x == pos_y
        % Reorder centroids based on robot order
        [Cv(:,pos_x), Cv_true(:,pos_x)] = dintegrate(vx,vy,p(i,:),ai(:,pos_x),pos_x);
    else
        disp('Mismatch in position found')
        [Cv(:,pos_x)] = p(i,:)';
        [Cv_true(:,pos_x)] = p(i,:)';
    end
end

% Laplacian: Shared edge length as weight
T = delaunayTriangulation(px,py);
ed = edges(T);
L = zeros(n,n);
for ind = 1:length(ed)
    r1 = ed(ind,1);
    r2 = ed(ind,2);
    points = intersect(c{order(r1)},c{order(r2)});
    if length(points)==2
        edge_len = pdist(v(points,:),'euclidean');
    else
        edge_len = 0;
    end
    L(r1,r2)= -edge_len;
    L(r2,r1)= -edge_len;
end
L = L + diag(-1*sum(L,2));
end