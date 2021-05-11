% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

function [Cvi,Cvi_true] = dintegrate(vx,vy,pi_t,ai_t,pos)
    global kappa K Fi a sens_info_flag;
    
    % Sample points
    cnt = 10000; % Sample count
    xp = min(vx) + rand(cnt,1)*(max(vx)-min(vx));
    yp = min(vy) + rand(cnt,1)*(max(vy)-min(vy));
    in = inpolygon(xp,yp,vx,vy); % subset points in voronoi region 
    xp = xp(in);
    yp = yp(in);
    
    % Integrals over voronoi region
    kq = kappa(xp,yp);
    if(sens_info_flag)
        phi_est = kq * a; % optimal configuration case when sensory info is known
    else 
        phi_est = kq * ai_t; % Case when sensory info is unknown and estimate is used
    end
    mv = sum(phi_est);
    cx = sum(xp.*phi_est)/mv;
    cy = sum(yp.*phi_est)/mv;
    phi = kq * a;
    cx_true = sum(xp.*phi)/sum(phi);    
    cy_true = sum(yp.*phi)/sum(phi);
    Cvi_true = [cx_true;cy_true]; 
    
    % Check for negative mass
    if mv < 0
       disp(strcat('Negative mass calculated: ',num2str(mv)));
       disp(pos)
       disp(ai_t')
    end
    
    % Check for calculated centroid in voronoi region
    if inpolygon(cx,cy,vx,vy) == 0
        disp('Centroid outside the region');
        disp(ai_t');
        Cvi = pi_t';
    else
        Cvi = [cx;cy];
    end
    
    % Update paramter Fi
    k1 = zeros(9,2);
    for i = 1:length(xp)
        q_pi = ([xp(i) yp(i)] - pi_t);
        k1 = k1 + kq(i,:)'*q_pi;
    end
    Fi(:,:,pos) = (1/mv)*(k1*K*k1');
   
end