%% -------------
% Switching Figure for the paper "A phyiscal model for exploring th dynamics of ultra-fast motions in mantis shrimp"
% Author : Nak-seung P. Hyun
% Date : 09/23/2020
% Output : Generating four subfigures 
%           (a) Lowerbound -Figure 1
%           (b) Upperbound -Figure 2
%           (c) Tendon force bound at equilibrium - Figure 3
%           (d) State space trajectory - Figure 4
%---------------
clear all;
% retreive parameters from setup
 
data_import = csvread('data_summary_PDP2.csv');
th1_i = data_import(:,1);
th2_i = data_import(:,2);
unit_rad2deg = 180/pi;

% Example number
exp_index = 15;

%% Generate Figue (a), (b), (d)

for contact_angle_flag = 1:2
    %% Initial condition setup
  
    if contact_angle_flag==1 
         contact_angle= th2_i(exp_index);
        
    elseif contact_angle_flag==2
         contact_angle= th2_i(exp_index)+ 7*pi/180;
    end
   
    init = [th1_i(exp_index); contact_angle;0;0];


    %% create theta_1 vector and dtheta_1 vectors for contour plot range ****
    th1_init = 110/unit_rad2deg;
    th1_final = 160/unit_rad2deg;
    th1_vec = linspace(th1_init,th1_final,1000);
    num_1 = length(th1_vec);
    dth1_vec = linspace(-1.5,1.5,1000);
    num_d1 = length(dth1_vec);

    [~,index_dth1_0] = min(abs(dth1_vec));
    % create mesh grid using th1 and dth1 vectors
    [X,Y] = meshgrid(th1_vec,dth1_vec);

    if contact_angle_flag==1
        data_phi3 = load('Figure_switching_data_phi3.mat');
        Feasible_lowerbound=data_phi3.Feasible_lowerbound;
        Feasible_upperbound=data_phi3.Feasible_upperbound;
%         F_equilibrium_3=data_phi3.F_equilibrium;
%         save('Figure_switching_data_phi3.mat','Feasible_lowerbound','Feasible_upperbound','F_equilibrium');
    elseif contact_angle_flag==2
        data_phi10 = load('Figure_switching_data_phi10.mat');
        Feasible_lowerbound=data_phi10.Feasible_lowerbound;
        Feasible_upperbound=data_phi10.Feasible_upperbound;
%         F_equilibrium_10=data_phi10.F_equilibrium;
%         save('Figure_switching_data_phi10.mat','Feasible_lowerbound','Feasible_upperbound','F_equilibrium');
    end

    %% ------- create Lower bound and upper bound for phi= 3 -------

    % 
    if contact_angle_flag ==1
        for ii=1:3
                        % choose levels to view on contour plot ****
            levels = linspace(0,17,17);
            n_lev = length(levels);

            if ii==1
                % Lower bound plot Figure X.(a)
                fignum = 1;
                bound_flag =1;
            elseif ii==2
                % Upper bound plot Figure X.(b)
                fignum = 2;
                bound_flag =2;
            elseif ii==3
                % Upper bound plot Figure X.(d) with trajectory
                fignum = 4;
                bound_flag =2;
            end
            figure(fignum)
            hold on;
            if bound_flag==1 % Plot lower bound for Figure (a)
                [C,h] =contourf(X*180/pi,Y*1e3,(Feasible_lowerbound),levels);                
            elseif bound_flag==2 % Plot upper bound for Figure(b, d)
                [C,h] =contourf(X*180/pi,Y*1e3,(Feasible_upperbound),levels);
            end
            
            xlabel('\boldmath$\theta_1$ (deg)','FontSize',16,'Interpreter','latex')
            ylabel('\boldmath$\dot{\theta_1}$ (rad/s)','FontSize',16,'Interpreter','latex')
            set(h,'LineColor','none')
            % create colormap
            map1 = brewermap(n_lev+2,'GnBu');
            map1 = map1(3:end,:);
            % create colorbar
            c1 = colorbar;
            colormap(map1)
            title(c1, '\boldmath$F_{sw}$ upper bound','FontSize',16,'Interpreter','latex')
            % title(c1, '\boldmath$F_{sw}$ upper
            % bound','FontSize',16,'Interpreter','latex')

            % Circle the initial condition
            plot(init(1)*180/pi,init(3),'ro','MarkerSize',6,'LineWidth',1.2);
            % Make a horizontal line for dth1=0
            plot([th1_vec(1)*unit_rad2deg,th1_vec(end)*unit_rad2deg],[0, 0],'b');


        %% plot various tendon velocity inputs
            if ii==3
            % colormap setup for trajectories
                map2 = brewermap(1,'RdPu');

            % calculate trajectory for various forces and plot
                filename = strcat('Simulation_data_experiment_14_k_s_63_k_h_750_testing_new_air_drag_1_k_arm_damp_4000_addedmass_ratio_18.mat');
                sim_test=load(filename);

                % retreive states
                th1 = sim_test.th1*unit_rad2deg;
                dth1 = sim_test.dth1*1e3;
                ind_switch=find(sim_test.t>=sim_test.T_switch,1);

                [~,ind_max_dth1]=max(sim_test.dth1);
                % find indices of switch, over-centering, v_max
            %     [ind_switch,ind_oc,ind_vm] = find_inds(t,states,params);

                % plot state space trajectory
                plot_color = map2(1,:); % get color for plot
                plot(th1(1:ind_max_dth1),dth1(1:ind_max_dth1),'Color',plot_color,'LineWidth',3.2)

                % plot indices of interest if they aren't at the first position
                if ind_switch > 1
                    p_s = plot(th1(ind_switch),dth1(ind_switch),'*','Color',plot_color,'MarkerSize',8,'LineWidth',1.2);
                end
                
                xlim([111,145]);
                ylim([-10,230]);
            end
        end
    end    
end

%% Generate Figue (c)
figure(3);
hold on;
xlimit = 160; % 160 degree
ylimit = 50;
% th2_sign_flip = 129.8;
for contact_angle_flag=1:2
    if contact_angle_flag==1
        data_phi3 = load('Figure_switching_data_phi3.mat');
        Feasible_lowerbound=data_phi3.Feasible_lowerbound;
        Feasible_upperbound=data_phi3.Feasible_upperbound;
        F_equilibrium=data_phi3.F_equilibrium;
        [~,index_equilibrium_limit] = min(abs(Feasible_upperbound(index_dth1_0,:)-F_equilibrium(index_dth1_0,:)));
        index_equilibrium_limit_3 =index_equilibrium_limit;
        index_dth1_0_3=index_dth1_0;
        F_equilibrium_3=F_equilibrium;
        colorcode_upper = 'b-';
        colorcode_force = 'r-';
        % Plot lower bound
        plot(th1_vec*unit_rad2deg,Feasible_lowerbound(index_dth1_0,:),'k'); 
        
    elseif contact_angle_flag==2
        data_phi10 = load('Figure_switching_data_phi10.mat');
        Feasible_lowerbound=data_phi10.Feasible_lowerbound;
        Feasible_upperbound=data_phi10.Feasible_upperbound;
        F_equilibrium=data_phi10.F_equilibrium;
        [~,index_equilibrium_limit] = min(abs(Feasible_upperbound(index_dth1_0,:)-F_equilibrium(index_dth1_0,:)));
        index_equilibrium_limit_10 =index_equilibrium_limit;
        index_dth1_0_10=index_dth1_0;
        F_equilibrium_10=F_equilibrium;
        colorcode_upper = 'b--';
        colorcode_force = 'r--';
    end
        % Plot upper bound
        plot(th1_vec*unit_rad2deg,Feasible_upperbound(index_dth1_0,:),colorcode_upper); % Upper bound plot
        plot(th1_vec*unit_rad2deg,F_equilibrium(index_dth1_0,:),colorcode_force);

end
plot(th1_vec(index_equilibrium_limit_3)*unit_rad2deg*[1, 1], [0, ylimit], 'k--');
plot(th1_vec(index_equilibrium_limit_10)*unit_rad2deg*[1, 1], [0, ylimit], 'k--');
% Equilibrium limit points (black dots)
plot(th1_vec(index_equilibrium_limit_3)*unit_rad2deg, F_equilibrium(index_dth1_0_3,index_equilibrium_limit_3),'ko','MarkerFaceColor','k','MarkerSize',8);
plot(th1_vec(index_equilibrium_limit_10)*unit_rad2deg, F_equilibrium_10(index_dth1_0_10,index_equilibrium_limit_10),'ko','MarkerFaceColor','k','MarkerSize',8);
% plot(th2_sign_flip*[1, 1], [0, ylimit], 'b--');
xlabel('\boldmath$\theta_1$ (deg)','FontSize',16,'Interpreter','latex')
ylabel('\boldmath$F$ (N)','FontSize',16,'Interpreter','latex')
xlim([th1_vec(1)*unit_rad2deg, xlimit])
ylim([-1, ylimit])
legend('lower bound', 'upperbound for \phi=3','equilibrium force for \phi=3','upperbound for \phi=10 ', 'equilibrium force for \phi=10')
%% ----------- functions -----------

function [th1_init] = calc_th1i(params)
% Function calculates the initial theta_1 value corresponding to the
% initial spring length. Takes in system parameters to calculate.

l_s0 = params.l_s0;
L_4 = params.L_4;
d_s = params.d_s;
alpha = params.alpha;

th1_init = acos((L_4^2 + d_s^2 - l_s0^2)/(2*L_4*d_s)) + alpha;
           
end