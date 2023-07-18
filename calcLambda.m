function lamda = calcLambda(q,qdot,qddot_prv,h,torque)

    %%% Value Assignments %%%

    i_w = h.i_w;
    m = h.m;
    r = h.r;
    b = h.b;
    I = h.I;
    theta = q(3);
    theta_dot = qdot(3);
    X_dot = qdot(1);
    Y_dot = qdot(2);
    tau1 = torque(1);
    tau2 = torque(2);
    tau3 = torque(3);
    tau4 = torque(4);
    stau = tau1+tau2+tau3+tau4;
    
    int_v = intVel(qddot_prv,qdot,q);

    d = (int_v(3)^2 + int_v(4)^2)/qdot(3);
    if isnan(d)
        d = b;
    end
    d_dot = (((int_v(1)*int_v(3) + int_v(2)*int_v(4)/sqrt(int_v(3)^2 + int_v(4)^2))*qdot(3))-(sqrt(int_v(3)^2 + int_v(4)^2)*qddot_prv(3)))/qdot(3)^2;

    K = d_dot

    %%% Left side matrix
    %First row
    inv_T1(1,1) = (i_w*(3*i_w + m*r^2))/(m*r^4 + 4*i_w*r^2);
    inv_T1(1,2) = (b*i_w)/(2*d*r^2);
    inv_T1(1,3) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(1,4) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(1,5) = -(b*i_w)/(2*d*r^2);
    inv_T1(1,6) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    
    %Second row
    
    inv_T1(2,1) = (b*i_w)/(2*d*r^2);
    inv_T1(2,2) = (d^2*m*r^2 + I*r^2+ 4*b^2*i_w)/(4*(d^2*r^2));
    inv_T1(2,3) = -(b*i_w)/(2*d*r^2);
    inv_T1(2,4) = -(b*i_w)/(2*d*r^2);
    inv_T1(2,5) = -(I*r^2 - d^2*m*r^2 + 4*b^2*i_w)/(4*(d^2*r^2));
    inv_T1(2,6) = (b*i_w)/(2*d*r^2);
    
    %Third row 
    
    inv_T1(3,1) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(3,2) = -(b*i_w)/(2*d*r^2);
    inv_T1(3,3) = (i_w*(3*i_w + m*r^2))/(m*r^4 + 4*i_w*r^2);
    inv_T1(3,4) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(3,5) = (b*i_w)/(2*d*r^2);
    inv_T1(3,6) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    
    % Fourth row
    
    inv_T1(4,1) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(4,2) = -(b*i_w)/(2*d*r^2);
    inv_T1(4,3) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(4,4) = (i_w*(3*i_w + m*r^2))/(m*r^4 + 4*i_w*r^2);
    inv_T1(4,5) = (b*i_w)/(2*d*r^2);
    inv_T1(4,6) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    
    % Fifth row
    
    inv_T1(5,1) = -(b*i_w)/(2*d*r^2);
    inv_T1(5,2) = -(I*r^2 - d^2*m*r^2 + 4*b^2*i_w)/(4*(d^2*r^2));
    inv_T1(5,3) = (b*i_w)/(2*d*r^2);
    inv_T1(5,4) = (b*i_w)/(2*d*r^2);
    inv_T1(5,5) = (d^2*m*r^2 + I*r^2 + 4*b^2*i_w)/(4*(d^2*r^2));
    inv_T1(5,6) = -(b*i_w)/(2*d*r^2);
    
    % Sixth row
    
    inv_T1(6,1) =  -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(6,2) = (b*i_w)/(2*d*r^2);
    inv_T1(6,3) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(6,4) = -(i_w*(i_w))/(m*r^4 + 4*i_w*r^2);
    inv_T1(6,5) = -(b*i_w)/(2*d*r^2);
    inv_T1(6,6) = (i_w*(3*i_w + m*r^2))/(m*r^4 + 4*i_w);

    %%% Right side matrix

    temp_1 = [ cos(theta)/m, sin(theta)/m, -b/I, r/i_w,     0,     0,     0;
    -sin(theta)/m, cos(theta)/m,  d/I,     0,     0,     0,     0;
     cos(theta)/m, sin(theta)/m,  b/I,     0, r/i_w,     0,     0;
     cos(theta)/m, sin(theta)/m,  b/I,     0,     0, r/i_w,     0;
    -sin(theta)/m, cos(theta)/m, -d/I,     0,     0,     0,     0;
     cos(theta)/m, sin(theta)/m, -b/I,     0,     0,     0, r/i_w];
    
    temp_2 = [cos(theta)*stau;
        sin(theta)*stau;
        b*(-tau4-tau1+tau2+tau3);
        tau1;
        tau2;
        tau3;
        tau4];
    
    
    temp_3 = [(-sin(theta)*theta_dot*X_dot) + (cos(theta)*theta_dot*Y_dot) ;
             (-cos(theta)*theta_dot*X_dot) + (-sin(theta)*theta_dot*Y_dot) + K*theta_dot;
             (-sin(theta)*theta_dot*X_dot) + (cos(theta)*theta_dot*Y_dot);
             (-sin(theta)*theta_dot*X_dot) + (cos(theta)*theta_dot*Y_dot);
             (-cos(theta)*theta_dot*X_dot) + (-sin(theta)*theta_dot*Y_dot) + K*theta_dot;
             (-sin(theta)*theta_dot*X_dot) + (cos(theta)*theta_dot*Y_dot)];

    T2= temp_1*temp_2 + temp_3


    %%% Lamda calculation

    lamda = inv_T1*T2;

end