function [drift_vec, actuate_vec, torque, torque_vec] = get_dynamics(z, Ft,params,M_q,C_q,G_q,B_q)


dth1=z(3);
dth3=z(4);
% [M_q, C_q, G_q, B_q] = get_Lagrangian(z,params);

e2=[0;1];
det_M_q = M_q(1,1)*M_q(2,2)-M_q(1,2)*M_q(2,1);
M_q_inv = 1/det_M_q*[M_q(2,2), -M_q(1,2);...
                     -M_q(2,1), M_q(1,1)];
drift_vec = M_q_inv*(-C_q*[dth1;dth3]-G_q);
actuate_vec = M_q_inv*B_q*Ft;

torque = -(e2'*(drift_vec+actuate_vec))/(e2'*(M_q_inv*e2));
torque_vec = M_q_inv*e2*torque;

end
