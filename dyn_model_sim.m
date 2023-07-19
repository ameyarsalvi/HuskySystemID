clc
clear



% let slip angle for each wheel be denoted by alpha_i

%alpha_i = invTan(V_y/V_x), 
% where V_y = x_projection_of_ICR * yaw_rate ;
% OR V_y can be integrated from IMU data
%V_x = phi_i_dot (can be calculated from wheel odometery)

% m = 81 Kgs

% I = I_cg + m*d^2 =  (1/12) m (a2 + b2) + m*d2 = 24.0253 + 81*0.0882 =
% 31.1695 kg*m^2

%I_w = 0.068062 kg*m^2

%% 

%%% Constants : 
h.m = 810; 
h.I = 311.695;
h.i_w = 0.68062;
h.r = 0.1650;
h.mu_lat = 2.5;
h.a = 0.6477;
h.b = 0.6858;
h.N_lf = 19.5*9.8;
h.N_lr = 22.5*9.8;
h.N_rf = 16*9.8;
h.N_rr = 23*9.8;


%%% State space

torque = [10;5;5;7];

q_prv = zeros(7,1);
qdot_prv = zeros(7,1);
qddot_prv = zeros(7,1);
theta = q_prv(3);
d=0;
qdot=[];
q =[];
qddot=[];
M_inv = [1/h.m,   0,   0,     0,     0,     0,     0;
  0, 1/h.m,   0,     0,     0,     0,     0;
  0,   0, 1/h.I,     0,     0,     0,     0;
  0,   0,   0, 1/h.i_w,     0,     0,     0;
  0,   0,   0,     0, 1/h.i_w,     0,     0;
  0,   0,   0,     0,     0, 1/h.i_w,     0;
  0,   0,   0,     0,     0,     0, 1/h.i_w];
 

% A_mat_prv = [cos(theta) sin(theta) -h.b h.r 0 0 0;
%              -sin(theta) cos(theta) d 0 0 0 0;
%              cos(theta) sin(theta) h.b 0 h.r 0 0;
%              cos(theta) sin(theta) h.b 0 0 h.r 0;
%              -sin(theta) cos(theta) -d 0 0 0 0;
%              cos(theta) sin(theta) -h.b 0 0 0 h.r];

dt =0.01;
t_end = 10;
i=1;

for t = 1:10

    [lambda, cr_mat, tau, R, M] = getMatrix(q_prv,qdot_prv,qddot_prv,h,torque,dt);

    p = tau - R + cr_mat
    qddot = M_inv*(tau - R + cr_mat);
    %qddot(i,:) = qddot; 
    qdot(i,:) = qdot_prv + (qddot).*dt;
    q(i,:) = q_prv + qdot(i,:)'.*dt;

    q_prv = q(end,:)';
    qdot_prv = qdot(end,:)';
    qddot_prv = qddot;
    %A_mat_prv = A_mat;
    i=i+1
end

%%