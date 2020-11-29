% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
function Cv = compute_centroid(px,py,ai)

global n;

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

    if pos_x == pos_y
        % Reorder centroids based on robot order
        [Cv(:,pos_x)] = dintegrate(vx,vy,p(i,:),ai(:,pos_x),pos_x);
    else
        disp('Mismatch in position found')
        [Cv(:,pos_x)] = p(i,:)';
    end
end