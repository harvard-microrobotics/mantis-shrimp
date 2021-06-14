%% setting up system parameters to pass

params.L_1 = 19; %20;
params.L_2 = 1.5;
params.L_3 = 18.66;
params.L_4 = 16.56;
params.L_5 = 5;

params.S_a = 121;
params.S_b = 300;

params.d_s = 18.31*1;
params.d_t = 23; %16.31*1;
params.l_s0 = 28*1;
%m_tip = .1*10^-3; % .1g
params.m_tip = 0.;

params.alpha = deg2rad(50.63);
params.beta = deg2rad(140.24);
params.angle_offset = deg2rad(0);

params.m_1 =  0.53; %0.306;
params.I_1 = 28.08*1e0; %16.04;
params.x_1 = 6.53; %4.42;
params.y_1 = 10.98; %11.03;
%COM_1 = [m_1;I_1;x_1;y_1];

params.m_2 = 0.115;
params.I_2 = 2.76; %2.77;
params.x_2 = 1.36; %1.38;%-3.83; % New coordinate
params.y_2 = 3.78;%3.38; %1.38; % New coordinate
%COM_2 = [m_2;I_2;x_2;y_2];

% hinge stiffnesses
params.kh_1 = 2*0.001465*10^3; % N-mm
params.kh_2 = 0.002302*10^3; % N-mm

% hinge initial position for th2
params.th2h_i = deg2rad(39.76);

% gravity constant
params.gravity = 9.807*10^-3;

% set contact offset
params.phi = deg2rad(3);

% Stiffness/ linear damping/ fluid density to the air  for th1 
params.k_s = 47.7;  
params.k_damp_linear = 65.30 ;
params.k_damp = 1; 

% Stiffness/ linear damping/ fluid density to the air  for th2 

params.kh_2 =  0.71;
params.k_arm_damp = 4.000; 
params.k_arm_damp_aero = 1;

%% Added mass set up
scale = 1e3;

params.rho = 1.225/scale^2;
params.epsilonn = 1e-4; % relaxation for th3 get close to pi

params.added_mass_flag = 1;
params.added_mass_geometry = 1; % 0 : flat plate , 1: eliptical
added_mass_aspect_ratio = 1/8;

rho_air_water = params.k_damp*params.rho; % params.k_damp =750 for water 1 for air

% arm blade length
e_2=[0;1];
R_m_beta = [cos(params.beta), sin(params.beta);...
         -sin(params.beta), cos(params.beta)];

added_mass_scale = 1/4;
added_mass_yaw_scale =1/16*1/8;

major_axis_arm = norm(params.L_3*e_2-R_m_beta*params.L_5*e_2);
minor_axis_arm = params.S_a/major_axis_arm;

added_mass_com = 1/2*(params.L_3*e_2-R_m_beta*params.L_5*e_2)-[params.x_2;params.y_2];

major_axis_body = params.L_1;
minor_axis_body = params.S_b/major_axis_body;

bottom_axis_body = 1/2;

added_mass = added_mass_scale*rho_air_water*pi*minor_axis_arm*major_axis_arm^2; %added_mass_scale*rho_water*pi*params.S_a*params.S_a/params.L_3;
added_mass_minor = added_mass_scale*rho_air_water*pi*minor_axis_arm*major_axis_arm^2*added_mass_aspect_ratio^2; %added_mass_scale*rho_water*pi*params.S_a*params.S_a/params.L_3;

added_mass_moment = added_mass_yaw_scale*rho_air_water*pi*minor_axis_arm*major_axis_arm^4*(1-added_mass_aspect_ratio^2)^2;%added_mass_yaw_scale*added_mass_scale*rho_water*pi*params.S_a*params.S_a*params.L_3^2;

added_mass_b = added_mass_scale*rho_air_water*pi*minor_axis_body*major_axis_body^2*(1+bottom_axis_body);%added_mass_scale*rho_water*pi*params.S_b*params.S_b/params.L_1;
added_mass_moment_b = added_mass_yaw_scale*rho_air_water*pi*minor_axis_body*major_axis_body^4*(1+bottom_axis_body^3);  %added_mass_yaw_scale*added_mass_scale*rho_water*pi*params.S_b*params.S_b*params.L_1;

params.scale_added_mass_length = 1.8;

% save it to params list
params.added_mass = added_mass;
params.added_mass_minor = added_mass_minor;
params.added_mass_moment = added_mass_moment;
params.added_mass_b = added_mass_b;
params.added_mass_moment_b = added_mass_moment_b;
