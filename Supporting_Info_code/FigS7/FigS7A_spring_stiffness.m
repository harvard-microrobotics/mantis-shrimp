%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for Figure S7 in "A Physical Model of Mantis Shrimp for
% Exploringthe Dynamics of Ultra-Fast Systems" PNAS
% Code Author : Nak-seung P. Hyun
% Date : April 11 2021
% Objective : Show the system tidentification of stiffness/ compare the input force from experimetn to simulation/ Analyze the different loading condition based on the hinge stiffness 
% Output : Fig.S7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
params_setup2_FigS7
fig_num=1;

%% Fig.S7A code generation

data=load('mode1_data_PDP.mat');

th1_eq = data.mode1.th_eq;
F_eq = data.mode1.F_eq;
th2_off = deg2rad(5); % offset angle during the Mode 1

valid_index = [1:12,16:size(th1_eq,1)];


for i=1:length(valid_index)
[~, ~, ~, B_q_after] = get_Lagrangian([th1_eq(valid_index(i),2),th2_off,0,0],params);
torque_after(i) = -B_q_after(1)*F_eq(valid_index(i),2);
end

Xt = [ones(length(valid_index),1)*(pi+params.alpha),th1_eq(valid_index,2)];
bt = Xt\[torque_after'];

% Xt = [th1_eq(:,2)-ones(size(th1_eq,1),1)*(pi+params.alpha) ; ...
%        th1_eq(:,1)-ones(size(th1_eq,1),1)*(pi+params.alpha)];
% bt = Xt\[torque_after';torque_before']

% Estimated stiffness
th_vec = linspace(min(min(Xt)),max(max(Xt)),100);
tau_vec = (th_vec)*bt(2)+bt(1)*(pi+params.alpha);

% Supplement Figure S7A 
figure(fig_num);
hold on;
% plot(th1_eq(valid_index,1),torque_before,'ro');
plot((pi+params.alpha)-th1_eq(valid_index,2),torque_after,'bo', 'MarkerFaceColor', 'b','MarkerSize',6);
plot((pi+params.alpha)-th_vec,tau_vec,'k');
xlabel('\boldmath$\Delta\theta_1$ (rad)','FontSize',16,'Interpreter','latex')
ylabel('Torque (Nmm)','FontSize',16,'Interpreter','latex')
xlim([1.4, 2])

%% Fig.S7B-C code generation
% Params setup
scale = 1e3;


L_1 = params.L_1;
L_2 = params.L_2;
L_3 = params.L_3;
L_4 = params.L_4;

% Load exp and sim data

force_data15_name = 'data_sysid_15_test' ;
force_data15 = load(force_data15_name);

file_name = 'data_set5.csv';
data_exp_15 = csvread(file_name);
t_exp5 = data_exp_15(:,1);
Ts = t_exp5(2)-t_exp5(1);
th1_exp5 = data_exp_15(:,2);
dth1_exp5 = data_exp_15(:,3)/scale;
th2_exp5 = data_exp_15(:,4);
dth2_exp5 = data_exp_15(:,5)/scale;
F_exp5 = force_data15.data.force_rs;

[vlp_num, vlp_den]=butter(3,1000*Ts);

dth1_exp5_filtered = filtfilt(vlp_num,vlp_den,dth1_exp5);
dth2_exp5_filtered = filtfilt(vlp_num,vlp_den,dth2_exp5);

ddth1_exp5 = filter([1 -1],1, dth1_exp5_filtered)/(Ts);
ddth2_exp5 = filter([1 -1],1, dth2_exp5_filtered)/(Ts);

% Find Phase 4 index
[~,ind_phase4_exp5] = min(th1_exp5);
list_min_phase4_exp5 = find(th1_exp5==th1_exp5(ind_phase4_exp5));
ind_phase4_exp5=list_min_phase4_exp5(end);

% Find the vtip max index
tend_exp5 = size(th1_exp5)-9270;
vtip_exp5 = calc_vtip(th2_exp5(1:tend_exp5),dth1_exp5(1:tend_exp5),dth2_exp5(1:tend_exp5),params);
[~,vmax_ind] = max(vtip_exp5);
t_exp5_vmax = t_exp5(vmax_ind);

% Plotting setup

vertical_phase4_line_th1 = linspace(min(th1_exp5),max(th1_exp5),100);

x_lim_lb = 2117;
x_lim_ub = 2123;

% Fig.S7B set up for Force plot of experiment plot
figure(fig_num+1)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale, linspace(0,18,100),'k--');
plot(t_exp5(1:vmax_ind)*scale,F_exp5(1:vmax_ind),'k');


% Fig.S7C set up for state space (th1, dth1) of experiment plot
figure(fig_num+2)
hold on
grid on

plot(th1_exp5(1:vmax_ind),dth1_exp5(1:vmax_ind),'k-','LineWidth',2)


for j=1:3
    if j==1 
        file_name_sim =strcat('Simulation_data_experiment_new_alpha_14_k_s_477_k_h_710_testing_new_air_drag_1_k_damp_linear_65300_k_arm_damp_4000_addedmass_ratio_18');
        colorcode = 'r';
    elseif j==2
        file_name_sim =strcat('Simulation_data_experiment_new_alpha_14_k_s_477_k_h_750_testing_new_air_drag_1_k_damp_linear_65300_k_arm_damp_4000_addedmass_ratio_18');
        colorcode = 'b-';
    elseif j==3
        file_name_sim =strcat('Simulation_data_experiment_new_alpha_14_k_s_477_k_h_650_testing_new_air_drag_1_k_damp_linear_65300_k_arm_damp_4000_addedmass_ratio_18');
        colorcode = 'g-';
    end
    data_sim_15 = load(file_name_sim);
    t_sim = data_sim_15.t;
    th1_sim = data_sim_15.th1;
    dth1_sim = data_sim_15.dth1;
    th2_sim = data_sim_15.th3;
    dth2_sim = data_sim_15.dth3;
    T_switch_sim = data_sim_15.T_switch
    
    F_sim = data_sim_15.F;

% Find the switching index from Mode 1 to Mode 2 in the simulation
if T_switch_sim~=0
    ind_switch_temp = find(t_sim>T_switch_sim,1);
    if isempty(ind_switch_temp)
        ind_switch = 1;
    else
    ind_switch=ind_switch_temp;
    
    end
else
    ind_switch=length(t_sim); 
end

% Find the start of Phase 4 in ths simulation
    [~,ind_phase4_sim] = min(th1_sim);

% Align the time to match the start of phase 4 with the experiment.
    tend_sim = size(t_sim)-570;
    t_offset(j) = t_exp5(ind_phase4_exp5)*scale -t_sim(ind_phase4_sim);
    
    vtip = calc_vtip(th2_sim(1:tend_sim),dth1_sim(1:tend_sim),dth2_sim(1:tend_sim),params);
    [~,vmax_ind_sim] = max(vtip);
    t_sim_vmax = t_sim(vmax_ind_sim);

    
% Plot Fig.S7B
    figure(fig_num+1)
    hold on;
    if j==1
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),F_sim(1:vmax_ind_sim),colorcode);
    end
    ylabel('\boldmath$F_t$ (N)','FontSize',16,'Interpreter','latex')
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')

    legend('Phase 4 start', 'Experiment','Simulation (k_s=47.7)')%, 'Sim (b_2=4)')% 

% Plot Fig.S7C    
    figure(fig_num+2)
    hold on
    plot(th1_sim(1:vmax_ind_sim),dth1_sim(1:vmax_ind_sim),colorcode,'LineWidth',2)
    xlabel('\boldmath$\theta_1$ (rad)','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath$\dot{\theta_1}$ (rad/s)','FontSize',16,'Interpreter','latex')
    legend('Experiments', 'Simulation (0.71)', 'Simulation (0.75)', 'Simulation (0.65)')

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
