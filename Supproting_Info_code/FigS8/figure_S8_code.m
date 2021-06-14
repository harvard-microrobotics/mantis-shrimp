%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for Figure S8 in "A Physical Model of Mantis Shrimp for
% Exploringthe Dynamics of Ultra-Fast Systems" PNAS
% Code Author : Nak-seung P. Hyun
% Date : April 11 2021
% Objective : Show the system tidentification of damping coefficient/ compare the state space trajectory of simiulation with multiple damping coefficients and the experiment data
% Output : Fig.S8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
params_setup2_FigS8
fig_num=1;

%%
%% Params setup
L_1 = params.L_1;
L_2 = params.L_2;
L_3 = params.L_3;
L_4 = params.L_4;

% params.alpha = deg2rad(27.9317)

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

% Find ending time
tend_exp5 = size(th2_exp5)-9270;

% Plotting setup
vertical_phase4_line_th2 = linspace(min(th2_exp5),max(th2_exp5),100);
vertical_phase4_line_th1 = linspace(min(th1_exp5),max(th1_exp5),100);

x_lim_lb = 2117;
x_lim_ub = 2123;


% Tendon length

% calculate tendon length
tendon_length_exp5 = zeros(length(t_exp5),1);
for i = 1:length(t_exp5)
    tendon_length_exp5(i) = calc_lt(th1_exp5(i),th2_exp5(i),params);
end
[vlp_num, vlp_den]=butter(3,30*Ts);
tendon_length_filtered = filtfilt(vlp_num,vlp_den,tendon_length_exp5);

dtendon_length =  filter([1 -1],1, tendon_length_filtered)/Ts;

% figure set up

% Fig.S8A
figure(fig_num)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale, vertical_phase4_line_th1,'k--');
plot(t_exp5(1:vmax_ind)*scale,th1_exp5(1:vmax_ind),'k');
xlim([x_lim_lb,x_lim_ub])

% Fig.S8B
figure(fig_num+1)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale, vertical_phase4_line_th2,'k--');
plot(t_exp5(1:vmax_ind)*scale,th2_exp5(1:vmax_ind),'k');
xlim([x_lim_lb,x_lim_ub])

% Fig.S8C
figure(fig_num+2)
hold on
grid on
plot(th1_exp5(1:vmax_ind),dth1_exp5(1:vmax_ind),'k-','LineWidth',2)

% Fig.S8D
figure(fig_num+3)
hold on
grid on
plot(th2_exp5(1:vmax_ind),dth2_exp5(1:vmax_ind),'k-','LineWidth',2)

% Fig.S8E
figure(fig_num+4)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale,linspace(min(tendon_length_exp5),max(tendon_length_exp5),100),'k--');

% Fig.S8F
figure(fig_num+5)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale, linspace(0,18,100),'k--');
plot(t_exp5(1:vmax_ind)*scale,F_exp5(1:vmax_ind),'k');



for j=1:3
    if j==1 
            file_name_sim =strcat('Simulation_data_experiment_new_alpha_14_k_s_477_k_h_710_testing_new_air_drag_1_k_damp_linear_65300_k_arm_damp_4000_addedmass_ratio_18');
        
        colorcode = 'r';
    elseif j==2
         file_name_sim =strcat('Simulation_data_experiment_new_alpha_14_k_s_477_k_h_710_testing_new_air_drag_1_k_damp_linear_65300_k_arm_damp_3000_addedmass_ratio_18');
        colorcode = 'b-';
    elseif j==3
         file_name_sim =strcat('Simulation_data_experiment_new_alpha_14_k_s_477_k_h_710_testing_new_air_drag_1_k_damp_linear_65300_k_arm_damp_2000_addedmass_ratio_18');
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
    
    
    % Find Phase 4 starting index
    [~,ind_phase4_sim] = min(th1_sim);
    
    tend_sim = size(t_sim)-570;
    t_offset(j) = t_exp5(ind_phase4_exp5)*scale -t_sim(ind_phase4_sim);
    
    % Find vtip max index
    vtip = calc_vtip(th2_sim(1:tend_sim),dth1_sim(1:tend_sim),dth2_sim(1:tend_sim),params);
    [~,vmax_ind_sim] = max(vtip);
    t_sim_vmax = t_sim(vmax_ind_sim);
    
    % Calculate tendon_length
    for i_lt = 1:length(t_sim)
    tendon_length_sim(i_lt) = calc_lt(th1_sim(i_lt),th2_sim(i_lt),params);
    end
    dtendon_length_sim(:,j) = filter([1 -1],1, tendon_length_sim)/Ts;
    
    
    % Fig.8A plot
    figure(fig_num)
    hold on;
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),th1_sim(1:vmax_ind_sim),colorcode);
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath$\theta_1$','FontSize',16,'Interpreter','latex')
    legend('Phase 4 start', 'Experiment','Sim (b_2=4)', 'Sim (b_2=3)', 'Sim (b_2=2)')
    
    % Fig.8B plot
    figure(fig_num+1)
    hold on;
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),th2_sim(1:vmax_ind_sim),colorcode);
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath$\theta_2$','FontSize',16,'Interpreter','latex')
    legend('Phase 4 start', 'Experiment','Sim (b_2=4)', 'Sim (b_2=3)', 'Sim (b_2=2)')

    % Fig.8C plot
    figure(fig_num+2)
    hold on
    plot(th1_sim(1:vmax_ind_sim),dth1_sim(1:vmax_ind_sim),colorcode,'LineWidth',2)
    xlabel('\boldmath$\theta_1$ (rad)','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath$\dot{\theta_1}$ (rad/s)','FontSize',16,'Interpreter','latex')
    legend('Experiments', 'Sim (b_2=4)', 'Sim (b_2=3)', 'Sim (b_2=2)')

    % Fig.8D plot
    figure(fig_num+3)
    hold on
    plot(th2_sim(1:vmax_ind_sim),dth2_sim(1:vmax_ind_sim),colorcode,'LineWidth',2)
    legend('Experiments', 'Sim (b_2=4)', 'Sim (b_2=3)', 'Sim (b_2=2)')
    xlabel('\boldmath$\theta_2$ (rad)','FontSize',16,'Interpreter','latex')
    ylabel('\boldmath$\dot{\theta_2}$ (rad/s)','FontSize',16,'Interpreter','latex')

    % Fig.8E plot
    figure(fig_num+4)
    hold on;
    plot(t_exp5(1:vmax_ind)*scale,tendon_length_exp5(1:vmax_ind),'k');
    if j==1
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),tendon_length_sim(1:vmax_ind_sim),colorcode);
    end
    ylabel('\boldmath$l_t$ (mm)','FontSize',16,'Interpreter','latex')
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')
    legend('Phase 4 start', 'Experiment', 'Simulation') %'Sim (b_2=4)')
    xlim([2065,2125])
    ylim([37.3,37.9])
    % Fig.8F 
    figure(fig_num+5)
    hold on;
    if j==1
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),F_sim(1:vmax_ind_sim),colorcode);
    end
    ylabel('\boldmath$F_t$ (N)','FontSize',16,'Interpreter','latex')
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')
    legend('Phase 4 start', 'Experiment','Simulation (k_s=47.7)')%, 'Sim (b_2=4)')% 
    xlim([2105,2125])
    ylim([2,17])
    clear tendon_length_sim dtendon_length_sim
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

