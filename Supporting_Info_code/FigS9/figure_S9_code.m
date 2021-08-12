%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for Figure S9 in "A Physical Model of Mantis Shrimp for
% Exploringthe Dynamics of Ultra-Fast Systems" PNAS
% Code Author : Nak-seung P. Hyun
% Date : April 11 2021
% Objective : Show the simulation results with different loading condition
% (Added mass) and compare with the experimental data.
% Output : Fig.S9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare with experiment and simulation
% clear all;
clear all;
close all;

params_setup2_FigS9;

% Parameter setup
L_1 = params.L_1;
L_2 = params.L_2;
L_3 = params.L_3;
L_4 = params.L_4;


fig_num=1;
scale = 1e3;
num_sim =1;
% Load exp and sim data

force_data15_name = 'data_water_3.mat' ;
force_data15 = load(force_data15_name);

file_name = '011620_19kfps_T_2X_W_2_1_2p3_-8_7_3_data.csv';

data_exp_15 = csvread(file_name);

truncate_end = size(data_exp_15(:,1))-1000;
t_exp5 = data_exp_15(1:truncate_end,1);
Ts = t_exp5(2)-t_exp5(1);
th1_exp5 = data_exp_15(1:truncate_end,2);
dth1_exp5 = data_exp_15(1:truncate_end,4)/scale;
th2_exp5 = data_exp_15(1:truncate_end,3);
dth2_exp5 = data_exp_15(1:truncate_end,5)/scale;
F_exp5 = force_data15.data.force_rs;
th2_init = th2_exp5(1);
% contact_angle = th2_init-deg2rad(2.75);
% params.contact_angle = contact_angle;

[vlp_num, vlp_den]=butter(3,1000*Ts);

dth1_exp5_filtered = filtfilt(vlp_num,vlp_den,dth1_exp5);
dth2_exp5_filtered = filtfilt(vlp_num,vlp_den,dth2_exp5);

ddth1_exp5 = filter([1 -1],1, dth1_exp5_filtered)/(Ts);
ddth2_exp5 = filter([1 -1],1, dth2_exp5_filtered)/(Ts);

% Tendon length

% calculate tendon length
tendon_length = zeros(length(t_exp5),1);
for i = 1:length(t_exp5)
    tendon_length(i) = calc_lt(th1_exp5(i),th2_exp5(i),params);
end
[vlp_num, vlp_den]=butter(3,30*Ts);
tendon_length_filtered = filtfilt(vlp_num,vlp_den,tendon_length);

dtendon_length =  filter([1 -1],1, tendon_length_filtered)/Ts;
% Phase 4 starting index
[~,ind_phase4_exp5] = min(th1_exp5);
list_min_phase4_exp5 = find(th1_exp5==th1_exp5(ind_phase4_exp5));
ind_phase4_exp5=list_min_phase4_exp5(end);

vertical_phase4_line_th2 = linspace(min(th2_exp5),max(th2_exp5),100);
vertical_phase4_line_th1 = linspace(min(th1_exp5),max(th1_exp5),100);

x_lim_lb = 2117;
x_lim_ub = 2123;

% Get vtip max index 
tend_exp5 = length(th1_exp5)-0;
vtip_exp5 = calc_vtip(th2_exp5(1:tend_exp5),dth1_exp5(1:tend_exp5),dth2_exp5(1:tend_exp5),params);
[~,vmax_ind] = max(vtip_exp5)
t_exp5_vmax = t_exp5(vmax_ind)
vmax_ind = vmax_ind;

% Fig.9A setup
figure(fig_num)
hold on
grid on
plot(th1_exp5(1:vmax_ind),dth1_exp5(1:vmax_ind),'k-','LineWidth',2)

% Fig.9B setup
figure(fig_num+1)
hold on
grid on
plot(th2_exp5(1:vmax_ind),dth2_exp5(1:vmax_ind),'k-','LineWidth',2)

% Fig.9C setup
figure(fig_num+2)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale,linspace(0,1,100),'k--');

% Fig.9D setup
figure(fig_num+3)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale,linspace(min(tendon_length),max(tendon_length),100),'k--');

% Plot simulation

for j=1:3
    if j==1
        file_name_sim =strcat('Simulation_new_alpha_14_c_angle_13_k_s_477_k_h_50_testing_new_water_nodrag_750_k_damp_linear_65300_k_arm_damp_4000_addedmass_ratio_10');
        
        colorcode = 'b-';
    
    elseif j==2
        file_name_sim =strcat('Simulation_new_alpha_14_c_angle_13_k_s_477_k_h_50_testing_new_water_drag_750_k_damp_linear_65300_k_arm_damp_4000_addedmass_ratio_10');
        colorcode = 'g';
    
    elseif j==3
        file_name_sim =strcat('Simulation_new_alpha_14_c_angle_13_k_s_477_k_h_50_testing_new_water_drag_750_k_damp_linear_65300_k_arm_damp_4000_addedmass_ratio_18');
        colorcode = 'r-';
    end
    data_sim_15 = load(file_name_sim);
    t_sim = data_sim_15.t;
    th1_sim = data_sim_15.th1;
    dth1_sim = data_sim_15.dth1;
    th2_sim = data_sim_15.th3;
    dth2_sim = data_sim_15.dth3;
    T_switch_sim = data_sim_15.T_switch
    
    F_sim = data_sim_15.F;
    
    % Find the start of Phase 4 in simulation
    [~,ind_phase4_sim] = min(th1_sim);
    
    tend_sim = size(t_sim)- 0; %570;
    t_offset(j) = t_exp5(ind_phase4_exp5)*scale -t_sim(ind_phase4_sim);
    
    % Find vtip max index
    vtip = calc_vtip(th2_sim(1:tend_sim),dth1_sim(1:tend_sim),dth2_sim(1:tend_sim),params);
    [~,vmax_ind_sim] = max(vtip);
    t_sim_vmax = t_sim(vmax_ind_sim);
    
    % Calculate the tendon length and the aerodynamic forces
    for i_lt = 1:length(t_sim)
    tendon_length_sim(i_lt) = calc_lt(th1_sim(i_lt),th2_sim(i_lt),params);
    [F_aero_lift, F_aero_drag]=get_aerodynamics_sysid([th1_sim(i_lt),th2_sim(i_lt), dth1_sim(i_lt), dth2_sim(i_lt)],params);
    F_aero_lift_list(:,i_lt)=F_aero_lift;
    F_aero_drag_list(:,i_lt)=F_aero_drag;
    F_aero_lift_norm(i_lt) = norm(F_aero_lift);
    F_aero_drag_norm(i_lt) = norm(F_aero_drag);
    end
    
    % Fig.9A plot
    figure(fig_num)
    hold on
    plot(th1_sim(1:vmax_ind_sim),dth1_sim(1:vmax_ind_sim),colorcode,'LineWidth',2)
    xlabel('\boldmath$\theta_1$ (rad)','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath$\dot{\theta_1}$ (rad/s)','FontSize',16,'Interpreter','latex')
    legend('Experiments', 'Sim (No added mass)', 'Sim (added mass \gamma=4.2)', 'Sim (added mass \gamma=13.8)')
    
    % Fig.9B plot
    figure(fig_num+1)
    hold on
    plot(th2_sim(1:vmax_ind_sim),dth2_sim(1:vmax_ind_sim),colorcode,'LineWidth',2)
    legend('Experiments', 'Sim (No added mass)', 'Sim (added mass \gamma=4.2)', 'Sim (added mass \gamma=13.8)')
    xlabel('\boldmath$\theta_2$ (rad)','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath$\dot{\theta_2}$ (rad/s)','FontSize',16,'Interpreter','latex')

        
    % Fig.9C plot
    figure(fig_num+2)
    hold on;
    if j==3
        plot(t_sim(1:vmax_ind_sim)+t_offset(j),F_aero_drag_norm(1:vmax_ind_sim),'r');
        plot(t_sim(1:vmax_ind_sim)+t_offset(j),F_aero_lift_norm(1:vmax_ind_sim),'k');
    end
    ylabel('\boldmath$F_{aero}$ (N)','FontSize',16,'Interpreter','latex')
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')
    legend('Phase 4 start','Drag', 'Lift')% 
   xlim([1855,1885])
    

    % Fig.9D plot
    figure(fig_num+3)
    hold on;
    if j==3
        plot(t_exp5(1:vmax_ind)*scale,tendon_length(1:vmax_ind),'k');
        plot(t_sim(1:vmax_ind_sim)+t_offset(j),tendon_length_sim(1:vmax_ind_sim),colorcode);
    end
    ylabel('\boldmath$l_t$ (mm)','FontSize',16,'Interpreter','latex')
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')
    xlim([1800,1900])
    ylim([35.2,36.2])
    legend('Phase 4 start', 'Experiment','Sim (added mass \gamma=13.8)')

end



function [v_tip] = calc_vtip(th3,dth1,dth3,params)
% Given the angle and it's angular velocity, calculate the velocity of the tip
%   takes in th3 dth1 dth3 and outputs tip velocity v_tip
L_1 = params.L_1;
L_3 = params.L_3;

eqn1 = L_1^2*dth1.^2 + L_3^2*(dth1-dth3).^2;
eqn2 = -2*L_1*L_3*dth1.*(dth1-dth3).*cos(th3);

v_tip = sqrt(eqn1+eqn2);

end


function [lt_out] = calc_lt(th1,th3,params)
% calculate the tendon length from the current angles
L_1 = params.L_1;
L_2 = params.L_2;
d_t = params.d_t;
beta = params.beta;
lt_out = sqrt(L_1^2+L_2^2+d_t^2-2*L_1*L_2*cos(th3+beta)-2*d_t*(L_1*cos(th1)-L_2*cos(th1-th3-beta)));
end
