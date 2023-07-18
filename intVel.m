function int_v = intVel(qddot_prv,qdot,q)

%This function returns the estimated values for internal velocity and accelerations
%which otherwise would be determined from the IMU readings

theta = q(3);
theta_dot = qdot(3);

qddot_prv(1)
qddot_prv(2)

acc_int = [-sin(theta)*theta_dot cos(theta)*theta_dot; -cos(theta)*theta_dot -sin(theta)*theta_dot]*[qddot_prv(1);qddot_prv(2)];

vel_int = [cos(theta) sin(theta); -sin(theta) cos(theta)]*[qdot(1);qdot(2)];

int_v = [acc_int vel_int]

end