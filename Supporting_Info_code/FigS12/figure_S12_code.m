%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for Figure S12 in "A Physical Model of Mantis Shrimp for
% Exploringthe Dynamics of Ultra-Fast Systems" PNAS
% Code Author : Nak-seung P. Hyun
% Date : April 11 2021
% Objective : Show the work done by the tendon and the net energy stored in
% the spring during the experiment.
% Output : Fig.S12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
params_setup2_FigS12
fig_num=1;


%% Fig.S12A-C code generation
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

% infinitesimal_work=F_exp5.*dtendon_length*Ts;
for i= 1:length(t_exp5)
    [M_q, C_q, G_q, B_q] = get_Lagrangian([th1_exp5(i),th2_exp5(i),dth1_exp5(i),dth2_exp5(i)],params);
     infinitesimal_work(i)=F_exp5(i)*[dth1_exp5(i), dth2_exp5(i)]*B_q*Ts;
     [drift_vec, actuate_vec, torque, torque_vec] = get_dynamics([th1_exp5(i),th2_exp5(i),dth1_exp5(i),dth2_exp5(i)], F_exp5(i),params,M_q,C_q,G_q,B_q);
     constraint=M_q*torque_vec;
     KE_exp5(i) = 1/2*[dth1_exp5(i), dth2_exp5(i)]*M_q*[dth1_exp5(i), dth2_exp5(i)]';
    if i==1
        work_by_tendon(i) = 0;
    else
%         infinitesimal_work(i) = [dth1_exp5(i), dth2_exp5(i)]*constraint*Ts;
        work_by_tendon(i) = work_by_tendon(i-1)+infinitesimal_work(i)*scale;
    end
end
% 

% Find Phase 4 index
[~,ind_phase4_exp5] = min(th1_exp5);
list_min_phase4_exp5 = find(th1_exp5==th1_exp5(ind_phase4_exp5));
ind_phase4_exp5=list_min_phase4_exp5(end);

% PE calculation of the experiment data
PE_exp5 = 1/2*params.k_s/1.0*(th1_exp5-(pi+params.alpha)).^2;

% plot setup for phase 4
vertical_phase4_line_th2 = linspace(min(th2_exp5),max(th2_exp5),100);
vertical_phase4_line_th1 = linspace(min(th1_exp5),max(th1_exp5),100);
vertical_phase4_line_Work = linspace(min(work_by_tendon),max(work_by_tendon),100);

% Find the vtip max index
tend_exp5 = size(th1_exp5)-9270;
vtip_exp5 = calc_vtip(th2_exp5(1:tend_exp5),dth1_exp5(1:tend_exp5),dth2_exp5(1:tend_exp5),params);
[~,vmax_ind] = max(vtip_exp5);
t_exp5_vmax = t_exp5(vmax_ind);

% Plotting setup

vertical_phase4_line_th1 = linspace(min(th1_exp5),max(th1_exp5),100);

x_lim_lb = 2117;
x_lim_ub = 2123;

% Fig.S12A setup and plot the experimental data
figure(fig_num)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale, vertical_phase4_line_Work,'k--');
plot(t_exp5(1:vmax_ind)*scale,work_by_tendon(1:vmax_ind),'r');
plot(t_exp5(1:vmax_ind)*scale, PE_exp5(1:vmax_ind)-PE_exp5(1),'b');

% Fig.S12B setup and plot the experimental data
figure(fig_num+1)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale, vertical_phase4_line_Work,'k--');
plot(t_exp5(1:vmax_ind)*scale,work_by_tendon(1:vmax_ind),'r');
plot(t_exp5(1:vmax_ind)*scale, PE_exp5(1:vmax_ind)-PE_exp5(1),'b');


% Fig.S7B set up for Force plot of experiment plot
figure(fig_num+2)
hold on;
plot(ones(100,1)*t_exp5(ind_phase4_exp5)*scale, linspace(0,18,100),'k--');
plot(t_exp5(1:vmax_ind)*scale,F_exp5(1:vmax_ind),'k');

for j=1:1
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
    T_switch_sim = data_sim_15.T_switch;
    
    F_sim = data_sim_15.F;
     PE_sim(:,j) = 1/2*params.k_s/1.0*(th1_sim-(pi+params.alpha)).^2;
   


    % Find the start of Phase 4 in ths simulation
    [~,ind_phase4_sim] = min(th1_sim);

    % Align the time to match the start of phase 4 with the experiment.
    tend_sim = size(t_sim)-570;
    t_offset(j) = t_exp5(ind_phase4_exp5)*scale -t_sim(ind_phase4_sim);
    
    vtip = calc_vtip(th2_sim(1:tend_sim),dth1_sim(1:tend_sim),dth2_sim(1:tend_sim),params);
    [~,vmax_ind_sim] = max(vtip);
    t_sim_vmax = t_sim(vmax_ind_sim);
    
    % Calculate the work done in the simulation.

    for i= 1:length(t_sim)
            [M_q, C_q, G_q, B_q] = get_Lagrangian([th1_sim(i),th2_sim(i),dth1_sim(i),dth2_sim(i)],params);
        infinitesimal_work_sim(i)=F_sim(i)*[dth1_sim(i), dth2_sim(i)]*B_q*Ts;
        if i==1
            work_by_tendon_sim(i,j) = 0;
        else
            work_by_tendon_sim(i,j) = work_by_tendon_sim(i-1)+infinitesimal_work_sim(i)*scale;
        end
    end
    
    
    % Plot Fig.S12A
    figure(fig_num)
    hold on;
    if j==1
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),work_by_tendon_sim(1:vmax_ind_sim,j),'r--');
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),PE_sim(1:vmax_ind_sim,j)-PE_sim(1,j),'b--');
    end
    ylabel('\boldmath$W$ (mJ)','FontSize',16,'Interpreter','latex')
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')
    legend('Phase 4 start', 'Tendon (experiment)','Spring (experiment)', 'Tendon (simulation)','Spring (simulation)')

    % Plot Fig.S12B
    figure(fig_num+1)
    hold on;
    if j==1
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),work_by_tendon_sim(1:vmax_ind_sim,j),'r--');
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),PE_sim(1:vmax_ind_sim,j)-PE_sim(1,j),'b--');
    end
    ylabel('\boldmath$W$ (mJ)','FontSize',16,'Interpreter','latex')
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')
    legend('Phase 4 start', 'Tendon (experiment)','Spring (experiment)', 'Tendon (simulation)','Spring (simulation)')
    xlim([2100,2125])
    ylim([10,38])
    
    % Plot Fig.S12C
    figure(fig_num+2)
    hold on;
    if j==1
    plot(t_sim(1:vmax_ind_sim)+t_offset(j),F_sim(1:vmax_ind_sim),colorcode);
    end
    ylabel('\boldmath$F_t$ (N)','FontSize',16,'Interpreter','latex')
    xlabel('Time (ms)','FontSize',16,'Interpreter','latex')

    legend('Phase 4 start', 'Experiment','Simulation (k_s=47.7)')%, 'Sim (b_2=4)')% 

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
