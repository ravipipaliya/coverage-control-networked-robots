function Cvi = dintegrate(vx,vy,pi_t,pos)
    global kappa K Li li Fi a ai w_t;
    cnt = 1000000;
    xp = min(vx) + rand(cnt,1)*(max(vx)-min(vx));
    yp = min(vy) + rand(cnt,1)*(max(vy)-min(vy));
    in = inpolygon(xp,yp,vx,vy);
    xp = xp(in);
    yp = yp(in);
%     ai_t = [100          0.1          0.1          0.1          0.1          0.1          0.1          0.1          100]';
    
    % Region integral
    kq = kappa(xp,yp);
    phi_est = kq * ai(:,pos);     
    mv = sum(phi_est);
    cx = sum(xp.*phi_est)/mv;
    cy = sum(yp.*phi_est)/mv;
%     disp(mv);
    if inpolygon(cx,cy,vx,vy) == 0 || mv < 0
        disp('Centroid outside the region');
        Cvi = pi_t';
    else
        Cvi = [cx;cy];
    end
    
    % True sensory info at robot location
%     size(tprev)
%     w_t = 0.0001;
    ki = kappa(pi_t(1),pi_t(2));
    phi_t = ki*a;
    Li(:,:,pos) = Li(:,:,pos) + w_t*(ki'*ki);
    li(:,pos) = li(:,pos) + w_t*(phi_t)*ki';
    
    k1 = zeros(9,2);
%     k2 = zeros(2,9);
    for i = 1:length(xp)
        q_pi = ([xp(i) yp(i)] - pi_t);
        k1 = k1 + kq(i,:)'*q_pi;
%         k2 = k2 + q_pi'*kq(i,:);
    end
    Fi(:,:,pos) = (1/mv)*(k1*K*k1');
    