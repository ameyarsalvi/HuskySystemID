function [lambda, cr_mat, tau, R, M] = getMatrix(q,qdot,qddot_prv,h,torque,A_mat_prv,dt)
    
    theta = q(3);

    %%% Matrix Definitions : 

    %Mass Matrix :
    M = [h.m 0 0 0 0 0 0;
        0 h.m 0 0 0 0 0;
        0 0 h.I 0 0 0 0;
        0 0 0 h.i_w 0 0 0;
        0 0 0 0 h.i_w 0 0;
        0 0 0 0 0 h.i_w 0;
        0 0 0 0 0 0 h.i_w];


    %%% Dynamic variables

    v_int = [cos(theta) sin(theta); -sin(theta) cos(theta)]*[qdot(1);qdot(2)];
    vx = v_int(1);
    vy = v_int(2);

    %Resistance Matrix
    Fr_X = -sin(theta)*h.mu_lat*(h.N_lf + h.N_lr + h.N_rf + h.N_rr)*sign(vy);
    Fr_Y = cos(theta)*h.mu_lat*(h.N_lf + h.N_lr + h.N_rf + h.N_rr)*sign(vy);
    M_r = h.a* (-h.mu_lat*(h.N_lr+h.N_rr) + h.mu_lat*(h.N_lf+h.N_rf))*sign(vy);
    
    R = [Fr_X; Fr_Y;M_r;0;0;0;0];
    
    %Torque Matrix
    tau_fl =torque(1);
    tau_fr = torque(2);
    tau_rr = torque(3);
    tau_rl = torque(4);

    tau_1 = cos(theta)*(tau_fl+tau_fr+tau_rr+tau_rl);
    tau_2 = sin(theta)*(tau_fl+tau_fr+tau_rr+tau_rl);
    tau_3 = h.b*(-tau_rl-tau_fl+tau_fr+tau_rr);
    
    tau =[tau_1;tau_2;tau_3;tau_fl;tau_fr;tau_rr;tau_rl];
    
    
    %Constraint Matrix
    

    d = (vx^2 + vy^2)/qdot(3)
    if isnan(d)
        d = h.b
    end
    
    A_mat = [cos(theta) sin(theta) -h.b h.r 0 0 0;
             -sin(theta) cos(theta) d 0 0 0 0;
             cos(theta) sin(theta) h.b 0 h.r 0 0;
             cos(theta) sin(theta) h.b 0 0 h.r 0;
             -sin(theta) cos(theta) -d 0 0 0 0;
             cos(theta) sin(theta) -h.b 0 0 0 h.r];

    A_t = [cos(theta) -sin(theta) cos(theta) cos(theta) -sin(theta) cos(theta);
    sin(theta) cos(theta) sin(theta) sin(theta) cos(theta) sin(theta);
    -h.b d h.b h.b -d -h.b;
    h.r 0 0 0 0 0;
    0 0 h.r 0 0 0;
    0 0 0 h.r 0 0;
    0 0 0 0 0 h.r];
    
    %A_dot = (A_mat -A_mat_prv)/dt;
    
    %k_dum = M\A_mat';
    %l_1 = A_mat*k_dum;
    %k_dum2 = M\tau ;
    %l_2 = (A_mat*k_dum2+ A_dot*qdot);
    
    %lambda1 = l_1\l_2;
    

    %cr_mat = lambda*A_mat';
    lambda = calcLambda(q,qdot,qddot_prv,h,torque)
    cr_mat = A_t * lambda;


end