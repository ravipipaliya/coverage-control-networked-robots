function Cvi = dintegrate(vx,vy,pi_t,ai_t,pos)
    global kappa K Fi;
    
    % Sample points in voronoi region
    cnt = 10000;
    xp = min(vx) + rand(cnt,1)*(max(vx)-min(vx));
    yp = min(vy) + rand(cnt,1)*(max(vy)-min(vy));
    in = inpolygon(xp,yp,vx,vy);
    xp = xp(in);
    yp = yp(in);
    
    % Region integral
    kq = kappa(xp,yp);
    phi_est = kq * ai_t;     
    mv = sum(phi_est);
    cx = sum(xp.*phi_est)/mv;
    cy = sum(yp.*phi_est)/mv;

    % Check
    if mv < 0
       disp(strcat('Negative mass calculated: ',num2str(mv)));
       disp(pos)
       disp(ai_t')
    end
    
    if inpolygon(cx,cy,vx,vy) == 0
        disp('Centroid outside the region');
        disp(ai_t');
        Cvi = pi_t';
    else
        Cvi = [cx;cy];
    end
    
    % Update Fi
    k1 = zeros(9,2);
    for i = 1:length(xp)
        q_pi = ([xp(i) yp(i)] - pi_t);
        k1 = k1 + kq(i,:)'*q_pi;
    end
    Fi(:,:,pos) = (1/mv)*(k1*K*k1');