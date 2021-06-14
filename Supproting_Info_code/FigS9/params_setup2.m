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

%arm blade length
