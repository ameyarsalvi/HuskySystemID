
clc
clear

%Define global symbolic variables

syms m I i_w theta b r d

M = [m 0 0 0 0 0 0;
    0 m 0 0 0 0 0;
    0 0 I 0 0 0 0;
    0 0 0 i_w 0 0 0;
    0 0 0 0 i_w 0 0;
    0 0 0 0 0 i_w 0;
    0 0 0 0 0 0 i_w];

A_t = [cos(theta) -sin(theta) cos(theta) cos(theta) -sin(theta) cos(theta);
    sin(theta) cos(theta) sin(theta) sin(theta) cos(theta) sin(theta);
    -b d b b -d -b;
    r 0 0 0 0 0;
    0 0 r 0 0 0;
    0 0 0 r 0 0;
    0 0 0 0 0 r];

A = [cos(theta) sin(theta) -b r 0 0 0;
    -sin(theta) cos(theta) d 0 0 0 0;
    cos(theta) sin(theta) b 0 r 0 0;
    cos(theta) sin(theta) b 0 0 r 0;
    -sin(theta) cos(theta) -d 0 0 0 0;
    cos(theta) sin(theta) -b 0 0 0 r];

%%

inv_M = inv(M);

Term_1 = A*inv_M*A_t;

inv(Term_1)

%%

% T1 is simplified version of term 1 as symb toolbox does not consider
% trigo rule

T1 = zeros(6,6);

%syms bs rs ds 


%First row
T1(1,1) = 1/m + b*b/I + r*r/i_w;
T1(1,2) = (-b*d)/I;
T1(1,3) = 1/m - b^2/I;
T1(1,4) = 1/m - b^2/I;
T1(1,5) = (b*d)/I;
T1(1,6) = 1/m + b^2/I;

%Second row

T1(2,1) = (-b*d)/I;
T1(2,2) = 1/m + d^2/I;
T1(2,3) = (b*d)/I;
T1(2,4) = (b*d)/I;
T1(2,5) = 1/m - d^2/I;
T1(2,6) = (-b*d)/I;

%Third row 

T1(3,1) = 1/m - b^2/I;
T1(3,2) = (b*d)/I;
T1(3,3) = 1/m + b^2/I + r^2/i_w;
T1(3,4) = 1/m + b^2/I;
T1(3,5) = (-b*d)/I;
T1(3,6) = 1/m - b^2/I;

% Fourth row

T1(4,1) = 1/m - b^2/I;
T1(4,2) = (b*d)/I;
T1(4,3) = 1/m + b^2/I;
T1(4,4) = 1/m + b^2/I + r^2/i_w;
T1(4,5) = (-b*d)/I;
T1(4,6) = 1/m - b^2/I;

% Fifth row

T1(5,1) = (b*d)/I;
T1(5,2) = 1/m - d^2/I;
T1(5,3) = (-b*d)/I;
T1(5,4) = (-b*d)/I;
T1(5,5) = 1/m + d^2/I;
T1(5,6) = (b*d)/I;

% Sixth row

T1(6,1) = 1/m + b^2/I;
T1(6,2) = (-b*d)/I;
T1(6,3) = 1/m - b^2/I;
T1(6,4) = 1/m - b^2/I;
T1(6,5) = (b*d)/I;
T1(6,6) = 1/m + b^2/I + r^2/i_w;


%%

inv_T1 = zeros(6,6);

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

%%
syms stau tau1 tau2 tau3 tau4 

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

rh_term = temp_1*temp_2;

syms theta_dot X_dot Y_dot K 

A_dot_qdot = [(-sin(theta)*theta_dot*X_dot) + (cos(theta)*theta_dot*Y_dot) ;
             (-cos(theta)*theta_dot*X_dot) + (-sin(theta)*theta_dot*Y_dot) + K*theta_dot;
             (-sin(theta)*theta_dot*X_dot) + (cos(theta)*theta_dot*Y_dot);
             (-sin(theta)*theta_dot*X_dot) + (cos(theta)*theta_dot*Y_dot);
             (-cos(theta)*theta_dot*X_dot) + (-sin(theta)*theta_dot*Y_dot) + K*theta_dot;
             (-sin(theta)*theta_dot*X_dot) + (cos(theta)*theta_dot*Y_dot)];

rh_term_final = rh_term  + A_dot_qdot;

%%

% lambda = 

lamda = inv_T1 * rh_term;

