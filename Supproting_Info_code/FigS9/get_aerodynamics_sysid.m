  function [F_aero_lift, F_aero_drag ] = get_aerodynamics_sysid(z,params)

%% system constants retrieval from parameter pass
global T_switch
L_1 = params.L_1;
L_2 = params.L_2;
L_3 = params.L_3;
L_5 = params.L_5;
%d_s = params.d_s;
d_t = params.d_t;

alpha = params.alpha;
beta = params.beta;

m_1 = params.m_1;
I_1 = params.I_1;
x_1 = params.x_1;
y_1 = params.y_1;

m_2 = params.m_2;
I_2 = params.I_2;
x_2 = params.x_2;
y_2 = params.y_2;

rho = params.rho;
S_a = params.S_a;
S_b = params.S_b;

%l_s0 = params.l_s0; % initial spring length
k_s = params.k_s; % spring constant
k_damp = params.k_damp;
k_damp_linear = params.k_damp_linear;

k_arm_damp =params.k_arm_damp;
k_arm_damp_aero=params.k_arm_damp_aero;
% hinge stiffnesses (added to G_q)
%kh_1 = params.kh_1; lumped with k_s
kh_2 = params.kh_2;
% alpha= deg2rad(27.9317);% params.alpha;
% gravity
g = params.gravity;

% retreive initial position
th2h_i = params.th2h_i;

% retrieve offset of contact
phi = params.phi;

%% turn states into variable names that correspond to what they are for reading equations easily
th1 = z(1);
th2 = z(2);
dth1 = z(3);
dth2 = z(4);

%% setup system of equations

%% Lagrangian derivation

R_th1 = [cos(th1-pi/2), -sin(th1-pi/2);...
         sin(th1-pi/2), cos(th1-pi/2)];
R_th3 = [cos(pi-th2), -sin(pi-th2);...
         sin(pi-th2), cos(pi-th2)];
R_beta = [cos(beta), -sin(beta);...
         sin(beta), cos(beta)];
R_m_beta = [cos(beta), sin(beta);...
         -sin(beta), cos(beta)];
R_th3_m_beta = [cos(pi-th2-beta), -sin(pi-th2-beta);...
         sin(pi-th2-beta), cos(pi-th2-beta)];
R_th1_p_th3 = [cos(th1-th2+pi/2), -sin(th1-th2+pi/2);...
         sin(th1-th2+pi/2), cos(th1-th2+pi/2)];
JJ = [0 , -1;...
      1, 0];
R_th1_JJ = [cos(th1), -sin(th1);...
         sin(th1), cos(th1)];
% R_th1_JJ_p_th3 = [cos(th1+pi-th2), -sin(th1+pi-th2);...
%          sin(th1+pi-th2), cos(th1+pi-th2)];
e_1=[1;0];
e_2=[0;1];

com1 = [x_1; y_1];
com2 = [x_2; y_2];
copa = 1/2*(L_3*e_2+R_m_beta*L_5*e_2);
blade_vec = R_th1_p_th3*(L_3*e_2-R_m_beta*L_5*e_2);
blade_vec_norm = blade_vec/norm(blade_vec);
blade_vec_norm_orth = JJ*blade_vec_norm;

J_aero = [L_1*R_th1*e_2+R_th1_p_th3*copa, R_th1_p_th3*(-copa)];
d_copa = J_aero*[dth1;dth2];
d_copa_effective = d_copa'*blade_vec_norm_orth;
% F_aero = -d_copa*abs(d_copa_effective)*k_arm_damp_aero*rho*S_a;

if -d_copa'*(blade_vec_norm_orth)>0
   d_copa_lift = -JJ*(-d_copa);
   if d_copa_lift'*blade_vec_norm_orth>0
      d_copa_lift=-d_copa_lift; 
   end
   d_copa_normal = blade_vec_norm_orth;
else
   d_copa_lift = JJ*(-d_copa);
   
   if d_copa_lift'*blade_vec_norm_orth>0
      d_copa_lift=-d_copa_lift; 
   end
   d_copa_normal = -blade_vec_norm_orth;
end

epsilonn = params.epsilonn;
if norm(d_copa)>epsilonn
    AoA = acos(d_copa'*blade_vec_norm/norm(d_copa));
    C_L0 = 1.8;
    C_D0 = 0.4;
    C_Dmax = 3.4;

    C_L = C_L0*sin(2*AoA);
    C_D = (C_D0+C_Dmax)/2 -(C_D0-C_Dmax)/2*cos(2*AoA);
else
    C_L =0;
    C_D =0;
end
F_aero_drag = 1/2*(-d_copa)*norm(d_copa)*k_arm_damp_aero*rho*S_a*C_D;
F_aero_lift = 1/2*(d_copa_lift)*norm(d_copa_lift)*k_arm_damp_aero*rho*S_a*abs(C_L);

F_aero = F_aero_drag+F_aero_lift;

% F_aero = 1/2*(d_copa_normal)*norm(d_copa)^2*k_arm_damp_aero*rho*S_a;
B_aero_th1 = [1;0]*(-1/2*L_1/2*S_b*rho*dth1*abs(dth1)*(L_1/2)^2);
B_aero_th2 = J_aero.'*(1/2*(-d_copa)*norm(d_copa)*rho*S_a*C_D+1/2*(d_copa_lift)*norm(d_copa_lift)*rho*S_a*abs(C_L));
tau_aero_body = -1/2*L_1/2*S_b*rho*dth1*abs(dth1)*(L_1/2)^2*k_damp;

end