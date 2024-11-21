clear;
clc;
close all;

%% User Input
% Modify the maximum specific growth rates of the intermediate strains
% (medium-producing [MED] and low-producing [LO]) values relative to that
% of the wild type (WT).

% MED maximum specific growth rate multiplier (valid values: [0.6,1]):
MED_multiplier = 0.7;

% LO maximum specific growth rate multipler (valid values: [0.6,1]):
LO_multiplier = 0.8;

% Show total biomass on plots
% Separate semibatch plot into biomass and betaC plots
% go up to 8 batches in semibatch - tune parameters - done;
% 50mL of shaking flask volume - done;
% actual batch at beginning of continuous
% cana try a 24hr batch period
% add a mutation rate as a function of O2 limitation - instead, did
% decrease of yield - done;
% induced low oxygen caused faster losses by a lot - done;

if (0.6 <= MED_multiplier && MED_multiplier <= 1) && (0.6 <= LO_multiplier && LO_multiplier <= 1)
    %% Setup
    mu_max_1 = 0.3466.*0.6;
    K_HI_G = 1;
    K_HI_O2 = 0.002./0.6;
    k_d_HI = 0.005;
    MR_HI = 0.00005;
    
    mu_max_2 = 0.3466.*LO_multiplier;
    K_LO_G = 1;
    K_LO_O2 = 0.002./LO_multiplier;
    k_d_LO = 0.005;
    MR_LO = 0.00005;
    
    mu_max_3 = 0.3466;
    K_NO_G = 1;
    K_NO_O2 = 0.002;
    k_d_NO = 0.005;
    
    Y_HI_G = 2;
    Y_MED_G = 2;
    Y_LO_G = 2;
    Y_NO_G = 2;
    
    Y_HI_O2 = 0.05;
    Y_MED_O2 = 0.05;
    Y_LO_O2 = 0.05;
    Y_NO_O2 = 0.05;

    Y_prod_G = 536.873/(40/6*180.156); % g-prod/g-glucose

    alpha_mult = 10;
    beta_mult = 1.2;
    
    alpha_HI_prod = 0.001.*alpha_mult;
    alpha_MED_prod = 0.00075.*alpha_mult;
    alpha_LO_prod = 0.0005.*alpha_mult;

    beta_HI_prod = 0.001.*beta_mult;
    beta_MED_prod = 0.00075.*beta_mult;
    beta_LO_prod = 0.0005.*beta_mult;
    
    mu_max_4 = 0.3466.*MED_multiplier;
    K_MED_G = 1;
    K_MED_O2 = 0.002./MED_multiplier;
    k_d_MED = 0.005;
    MR_MED = 0.00005;
    
    p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
            mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
            alpha_HI_prod,alpha_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,alpha_MED_prod, ...
            beta_HI_prod,beta_MED_prod,beta_LO_prod,Y_prod_G];
    
    param_names = "mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO,mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2,alpha_HI_prod,alpha_LO_prod," + ...
        "mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,alpha_MED_prod,beta_HI_prod,beta_MED_prod,beta_LO_prod,Y_prod_G";
    param_names = split(param_names,',');
    
    disp(table((1:1:length(p))',p',param_names,'VariableNames',{'Number','Value','Name'}))
    
    % Environmental Conditions
    dV = 0.06; % L/s
    V = 1.5;% L
    C_G_in = 20; % g/L
    S_O2 = 0.007;
    k_L_a = 10;
    % k_L_a = 2;
    % amp = 0.075;
    amp = 0;
    freq = 50;
    y_shift = 0.3;
    env_cond = [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift];
    
    %% Run Model
    % Task 1
    % MED_multipler = 0.8
    % LO_multiplier = 0.9
    fig1_a = figure;
    fig1_a.Color = [1,1,1];
    
    cultivation_time = 300; % hr
    t_span = [0,cultivation_time];
    
    % base dilution rate (0.04 1/hr)
    ODE_1 = ode;
    ODE_1.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
    ODE_1.InitialTime = 0;
    ODE_1.InitialValue = [1.25,0.25,0.05,20,0.007,0,0];
    ODE_1.Parameters = p;
    ODE_1.Solver = 'ode15s';
    
    sol = solve(ODE_1,t_span(1),t_span(2));
    t_res_1_a = sol.Time;
    y_res_1_a = sol.Solution;
    
    % subplot(1,2,1);
    hold on;
    % yyaxis left;
    plot(t_res_1_a,y_res_1_a(1,:),'r-','LineWidth',1.5);
    plot(t_res_1_a,y_res_1_a(2,:),'k-','LineWidth',1.5);
    plot(t_res_1_a,y_res_1_a(3,:),'b-','LineWidth',1.5);
    plot(t_res_1_a,sum(y_res_1_a(1:3,:),1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
    yyaxis right;
    plot(t_res_1_a,y_res_1_a(5,:),'r--','LineWidth',1.5);
    axes = fig1_a.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err
    
        end
    end
    % yyaxis right;
    % plot(t_res_1_a,y_res_1_a(4,:),'b--','LineWidth',1.5);
    % ylim([0,0.15.*1000]);
    % legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass"],'Location','best');
    xlabel("Time (hr)"); % ylabel("\beta-Carotene Concentration (mg/L)");
    title("1_a");
    hold off;
    
    % subplot(3,2,3);
    % hold on;
    % yyaxis left;
    % plot(t_res,y_res(4,:),'b-','LineWidth',1.5);
    % ylabel("Glucose");
    % axes = fig1.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % plot(t_res,y_res(6,:),'k-','LineWidth',1.5);
    % % plot(t_res,y_res(7,:),'r-','LineWidth',1.5);
    % ylim([0,0.5]);
    % legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (mg/L)"],'Location','best');
    % xlabel("Time (hr)"); ylabel("\beta-Carotene");
    % hold off;
    % 
    % subplot(3,2,5);
    % hold on;
    % plot(t_res,y_res(5,:),'b-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
    % ylim([0,Inf]);
    % hold off;
    
    % higher dilution rate (0.08 1/hr)
    fig1_d = figure();
    fig1_d.Color = [1,1,1];

    env_cond(1) = 0.12; % 1/hr
    ODE_2 = ode;
    ODE_2.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
    ODE_2.InitialTime = 0;
    ODE_2.InitialValue = [1.25,0.25,0.05,20,0.007,0,0];
    ODE_2.Parameters = p;
    ODE_2.Solver = 'ode15s';
    env_cond(1) = 0.06; % 1/hr
    
    sol = solve(ODE_2,t_span(1),t_span(2));
    t_res_1_d = sol.Time;
    y_res_1_d = sol.Solution;
    
    % subplot(1,2,2);
    hold on;
    yyaxis left;
    plot(t_res_1_d,y_res_1_d(1,:),'r-','LineWidth',1.5);
    plot(t_res_1_d,y_res_1_d(2,:),'k-','LineWidth',1.5);
    plot(t_res_1_d,y_res_1_d(3,:),'b-','LineWidth',1.5);
    plot(t_res_1_d,sum(y_res_1_d(1:3,:),1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
    axes = fig1_d.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err
    
        end
    end
    % yyaxis right;
    % plot(t_res_1_d,y_res_1_d(6,:).*1000,'b--','LineWidth',1.5);
    % ylim([0,0.15.*1000]);
    % legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass"],'Location','best');
    xlabel("Time (hr)"); % ylabel("\beta-Carotene Concentration (mg/L)");
    title("Dilution Rate: 0.1 1/{hr} 1_d");
    hold off;
    
    % subplot(3,2,4);
    % hold on;
    % yyaxis left;
    % plot(t_res_2,y_res_2(4,:),'b-','LineWidth',1.5);
    % ylabel("Glucose"); ylim([0,Inf]);
    % axes = fig1.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % plot(t_res_2,y_res_2(6,:),'k-','LineWidth',1.5);
    % % plot(t_res_2,y_res_2(7,:),'r-','LineWidth',1.5);
    % ylim([0,0.5]);
    % legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (mg/L)"],'Location','best');
    % xlabel("Time (hr)"); ylabel("\beta-Carotene");
    % hold off;
    % 
    % subplot(3,2,6);
    % hold on;
    % plot(t_res_2,y_res_2(5,:),'b-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
    % ylim([0,Inf]);
    % hold off;
    
    axes = fig1_a.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err
    
        end
    end

    axes = fig1_d.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err
    
        end
    end

    disp(table(y_res_1_a(7,end).*1000,y_res_1_d(7,end).*1000,'VariableNames',{'0.05 1/hr','0.1 1/hr'}))
    
    % Task 2
    fig2a = figure;
    fig2a.Color = [1,1,1];
    
    cultivation_time = 300; % hr
    t_span = [0,cultivation_time];
    
    % 2 different strains
    ODE_3 = ode;
    ODE_3.ODEFcn = @(t,y,p) ODESys_cont_2(t,y,p,env_cond);
    ODE_3.InitialTime = 0;
    ODE_3.InitialValue = [1.5,0,0.05,20,0.007,0,0];
    ODE_3.Parameters = p;
    ODE_3.Solver = 'ode15s';
    
    sol = solve(ODE_3,t_span(1),t_span(2));
    t_res_2_a = sol.Time;
    y_res_2_a = sol.Solution;
    
    hold on;
    % yyaxis left;
    plot(t_res_2_a,y_res_2_a(1,:),'r-','LineWidth',1.5);
    plot(t_res_2_a,y_res_2_a(3,:),'b-','LineWidth',1.5);
    plot(t_res_2_a,sum(y_res_2_a(1:3,:),1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
    axes = fig2a.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    yyaxis right;
    plot(t_res_2_a,y_res_2_a(6,:).*1000,'r-.','LineWidth',3);
    ylim([0,0.15.*1000]);
    legend(["High-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    % legend(["High-Producer","Non-Producer","Total Biomass"],'Location','best');
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)");
    title("2 Distinct Strains 2_a");
    hold off;
    axes = fig2a.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    
    % subplot(3,3,4);
    % hold on;
    % yyaxis left;
    % plot(t_res,y_res(4,:),'b-','LineWidth',1.5);
    % ylabel("Glucose");
    % axes = fig2.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % plot(t_res,y_res(6,:),'k-','LineWidth',1.5);
    % % plot(t_res,y_res(7,:),'r-','LineWidth',1.5);
    % ylim([0,0.5]);
    % legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (mg/L)"],'Location','best');
    % xlabel("Time (hr)"); ylabel("\beta-Carotene");
    % hold off;
    % 
    % subplot(2,3,4);
    % hold on;
    % plot(t_res,y_res(5,:),'b-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
    % ylim([0,Inf]);
    % hold off;
    
    % 3 different strains
    fig2_bc = figure;
    fig2_bc.Color = [1,1,1];

    ODE_4 = ode;
    ODE_4.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
    ODE_4.InitialTime = 0;
    ODE_4.InitialValue = [1.25,0.25,0.05,20,0.007,0,0];
    ODE_4.Parameters = p;
    ODE_4.Solver = 'ode15s';
    
    sol = solve(ODE_4,t_span(1),t_span(2));
    t_res_2_bc = sol.Time;
    y_res_2_bc = sol.Solution;
    
    subplot(1,2,1);
    hold on;
    % yyaxis left;
    plot(t_res_2_bc,y_res_2_bc(1,:),'r-','LineWidth',1.5);
    plot(t_res_2_bc,y_res_2_bc(2,:),'k-','LineWidth',1.5);
    plot(t_res_2_bc,y_res_2_bc(3,:),'b-','LineWidth',1.5);
    plot(t_res_2_bc,sum(y_res_2_bc(1:3,:),1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
    axes = fig2_bc.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    yyaxis right;
    plot(t_res_2_bc,y_res_2_bc(6,:).*1000,'r-.','LineWidth',3);
    ylim([0,0.15.*1000]);
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)")
    title("3 Distinct Strains");
    hold off;
    axes = fig2_bc.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    
    % subplot(3,3,5);
    % hold on;
    % yyaxis left;
    % plot(t_res_2,y_res_2(4,:),'b-','LineWidth',1.5);
    % ylabel("Glucose"); ylim([0,Inf]);
    % axes = fig2.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % plot(t_res_2,y_res_2(6,:),'k-','LineWidth',1.5);
    % % plot(t_res_2,y_res_2(7,:),'r-','LineWidth',1.5);
    % ylim([0,0.5]);
    % legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (mg/L)"],'Location','best');
    % xlabel("Time (hr)"); ylabel("\beta-Carotene");
    % hold off;
    % 
    % subplot(2,3,5);
    % hold on;
    % plot(t_res_2,y_res_2(5,:),'b-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
    % ylim([0,Inf]);
    % hold off;
    
    % 4 different strains
    ODE_5 = ode;
    ODE_5.ODEFcn = @(t,y,p) ODESys_cont_4(t,y,p,env_cond);
    ODE_5.InitialTime = 0;
    ODE_5.InitialValue = [1,0.25,0.05,20,0.007,0,0,0.25];
    ODE_5.Parameters = p;
    ODE_5.Solver = 'ode15s';
    
    sol = solve(ODE_5,t_span(1),t_span(2));
    t_res_de = sol.Time;
    y_res_de = sol.Solution;
    
    subplot(1,2,2);
    hold on;
    % yyaxis left;
    plot(t_res_de,y_res_de(1,:),'r-','LineWidth',1.5);
    plot(t_res_de,y_res_de(8,:),'g-','LineWidth',1.5);
    plot(t_res_de,y_res_de(2,:),'k-','LineWidth',1.5);
    plot(t_res_de,y_res_de(3,:),'b-','LineWidth',1.5);
    plot(t_res_de,sum([y_res_de(1:3,:);y_res_de(8,:)],1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
    axes = fig2_bc.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    yyaxis right;
    plot(t_res_de,y_res_de(6,:).*1000,'r-.','LineWidth',3);
    ylim([0,0.15.*1000]);
    legend(["High-Producer","Medium-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)");
    title("4 Distinct Strains");
    hold off;
    axes = fig2_bc.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    
    % subplot(3,3,6);
    % hold on;
    % yyaxis left;
    % plot(t_res_3,y_res_3(4,:),'b-','LineWidth',1.5);
    % ylabel("Glucose"); ylim([0,Inf]);
    % axes = fig2.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % plot(t_res_3,y_res_3(6,:),'k-','LineWidth',1.5);
    % % plot(t_res_3,y_res_3(7,:),'r-','LineWidth',1.5);
    % ylim([0,0.5]);
    % legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (mg/L)"],'Location','best');
    % xlabel("Time (hr)"); ylabel("\beta-Carotene");
    % hold off;
    % 
    % subplot(2,3,6);
    % hold on;
    % plot(t_res_3,y_res_3(5,:),'b-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
    % ylim([0,Inf]);
    % hold off;
    
    axes = fig2a.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err
    
        end
    end

    % 3 different strains - no mutation
    fig1_c = figure;
    fig1_c.Color = [1,1,1];

    ODE_8 = ode;
    ODE_8.ODEFcn = @(t,y,p) ODESys_cont_3_no_mut(t,y,p,env_cond);
    ODE_8.InitialTime = 0;
    ODE_8.InitialValue = [1.25,0.25,0.05,20,0.007,0,0];
    ODE_8.Parameters = p;
    ODE_8.Solver = 'ode15s';
    
    sol = solve(ODE_8,t_span(1),t_span(2));
    t_res_4 = sol.Time;
    y_res_4 = sol.Solution;
    
    hold on;
    % yyaxis left;
    plot(t_res_4,y_res_4(1,:),'r-','LineWidth',1.5);
    plot(t_res_4,y_res_4(2,:),'k-','LineWidth',1.5);
    plot(t_res_4,y_res_4(3,:),'b-','LineWidth',1.5);
    plot(t_res_4,sum(y_res_4(1:3,:),1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
    axes = fig1_c.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    yyaxis right;
    plot(t_res_4,y_res_4(6,:).*1000,'r--','LineWidth',1.5);
    ylim([0,0.15.*1000]);
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)")
    title("3 Distinct Strains - No Mutation 1_c");
    hold off;

    axes = fig1_c.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    
    % subplot(3,3,5);
    % hold on;
    % yyaxis left;
    % plot(t_res_2,y_res_2(4,:),'b-','LineWidth',1.5);
    % ylabel("Glucose"); ylim([0,Inf]);
    % axes = fig2.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % plot(t_res_2,y_res_2(6,:),'k-','LineWidth',1.5);
    % % plot(t_res_2,y_res_2(7,:),'r-','LineWidth',1.5);
    % ylim([0,0.5]);
    % legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (mg/L)"],'Location','best');
    % xlabel("Time (hr)"); ylabel("\beta-Carotene");
    % hold off;
    % 
    % subplot(2,3,5);
    % hold on;
    % plot(t_res_4,y_res_4(5,:),'b-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
    % ylim([0,Inf]);
    % hold off;
    
    disp(table(y_res_2_a(7,end).*1000,y_res_2_bc(7,end).*1000,y_res_de(7,end).*1000,y_res_4(7,end).*1000,'VariableNames',{'2','3','4','3 - no mut'}))
    
    % Task 3
    % fig3 = figure;
    % fig3.Color = [1,1,1];
    % ODE_6 = ode;
    % ODE_6.ODEFcn = @(t,y,p) simple_mu_exper(t,y,p);
    % ODE_6.InitialTime = 0;
    % ODE_6.InitialValue = [1,0.01,20];
    % params = [[1,2].*0.08,[1,1,0.5,0.5]];
    % ODE_6.Parameters = params;
    % ODE_6.Solver = 'ode15s';
    % 
    % t_span = [0,300];
    % sol = solve(ODE_6,t_span(1),t_span(2));
    % t_res = sol.Time;
    % y_res = sol.Solution;
    % 
    % subplot(1,3,1);
    % hold on;
    % plot(t_res,ones(size(t_res)).*params(1),'b-','LineWidth',1.5);
    % plot(t_res,ones(size(t_res)).*params(2),'k-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("Specific Maximum Growth Rate (1/hr)");
    % legend(["ES","WT"],'Location','northeast');
    % ylim([0,0.3]);
    % hold off;
    % 
    % subplot(1,3,2);
    % hold on;
    % plot(t_res,params(1).*y_res(1,:).*(y_res(3,:)./(p(3)+y_res(3,:))),'b-','LineWidth',1.5);
    % plot(t_res,params(2).*y_res(2,:).*(y_res(3,:)./(p(4)+y_res(3,:))),'k-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("Growth Rate (g/(L*hr))");
    % legend(["ES","WT"],'Location','best');
    % hold off;
    % 
    % subplot(1,3,3);
    % hold on;
    % plot(t_res,y_res(1,:),'b-','LineWidth',1.5);
    % plot(t_res,y_res(2,:),'k-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
    % legend(["ES","WT"],'Location','best');
    % hold off;
    % 
    % axes = fig3.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % 
    % % closer mu
    % fig4 = figure;
    % fig4.Color = [1,1,1];
    % ODE_6 = ode;
    % ODE_6.ODEFcn = @(t,y,p) simple_mu_exper(t,y,p);
    % ODE_6.InitialTime = 0;
    % ODE_6.InitialValue = [1,0.01,20];
    % params = [[1.5,2].*0.08,[1,1,0.5,0.5]];
    % ODE_6.Parameters = params;
    % ODE_6.Solver = 'ode15s';
    % 
    % sol = solve(ODE_6,t_span(1),t_span(2));
    % t_res = sol.Time;
    % y_res = sol.Solution;
    % 
    % subplot(1,3,1);
    % hold on;
    % plot(t_res,ones(size(t_res)).*params(1),'b-','LineWidth',1.5);
    % plot(t_res,ones(size(t_res)).*params(2),'k-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("Specific Maximum Growth Rate (1/hr)");
    % legend(["ES","WT"],'Location','northeast');
    % ylim([0,0.3]);
    % hold off;
    % 
    % subplot(1,3,2);
    % hold on;
    % plot(t_res,params(1).*y_res(1,:).*(y_res(3,:)./(p(3)+y_res(3,:))),'b-','LineWidth',1.5);
    % plot(t_res,params(2).*y_res(2,:).*(y_res(3,:)./(p(4)+y_res(3,:))),'k-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("Growth Rate (g/(L*hr))");
    % legend(["ES","WT"],'Location','best');
    % ylim([0,3.5]);
    % hold off;
    % 
    % subplot(1,3,3);
    % hold on;
    % plot(t_res,y_res(1,:),'b-','LineWidth',1.5);
    % plot(t_res,y_res(2,:),'k-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
    % legend(["ES","WT"],'Location','best');
    % hold off;
    % 
    % axes = fig4.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end

    % Task 4
    batch_ct = 1;
    t_res_1_b = [];
    y_res_1_b = cell(6,1);
    tot_prod_mass = 0;
    V = 0.05; % L
    batch_time_length = 24; % hr
    t_span_1 = [0,batch_time_length];
    ODE_7 = ode;
    ODE_7.ODEFcn = @(t,y,p) ODESys_batch_3(t,y,p,env_cond);
    ODE_7.InitialTime = 0;
    ODE_7.InitialValue = [0.1,0.05,0,20,0.007,0];
    ODE_7.Parameters = p;
    ODE_7.Solver = 'ode15s';
    % evt = odeEvent(EventFcn=@(t,y,p) EvtFcn(t,y,p),Response="stop");
    % ODE_7.EventDefinition = evt;
    while (batch_ct < 9)
        sol = solve(ODE_7,t_span_1(1),t_span_1(2));

        t_res_1_b = [t_res_1_b;sol.Time']; %#ok<*AGROW>
        for l=1:1:size(sol.Solution,1)
            y_res_1_b{l} = [y_res_1_b{l},sol.Solution(l,:)]; %#ok<*SAGROW>
        end

        t_span_1(1) = sol.Time(end);
        ODE_7.InitialTime = t_span_1(1);
        t_span_1(2) = t_span_1(1)+batch_time_length;

        ODE_7.InitialValue = [y_res_1_b{1}(end).*0.01,y_res_1_b{2}(end).*0.01,y_res_1_b{3}(end).*0.01,20,0.007,0];
        % ODE_7.InitialValue = [0.125,0.025,0.005,20,0.007,0];

        tot_prod_mass(end+1) = y_res_1_b{6}(end).*V;

        batch_ct = batch_ct + 1;
    end

    tot_prod_mass_cum = zeros(size(tot_prod_mass));
    for k=1:1:length(tot_prod_mass)
        tot_prod_mass_cum(k) = sum(tot_prod_mass(1:k));
    end

    fig1_b = figure;
    fig1_b.Color = [1,1,1];

    % subplot(2,2,1);
    hold on;
    % yyaxis left;
    plot(t_res_1_b,y_res_1_b{1},'r-','LineWidth',4);
    plot(t_res_1_b,y_res_1_b{2},'k-','LineWidth',1.5);
    plot(t_res_1_b,y_res_1_b{3},'b-','LineWidth',1.5);
    plot(t_res_1_b,sum([y_res_1_b{1};y_res_1_b{2};y_res_1_b{3}],1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    axes = fig1_b.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    % yyaxis right;
    % plot(t_res_1_b,y_res_1_b{4},'b--','LineWidth',1.5);
    % ylim([0,0.15.*1000])
    % legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"]);
    xlim([0,200]);
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass"]);
    xlabel("Time (hr)"); % ylabel("\beta-Carotene Concentration (mg/L)");
    title("Repeating Batch Culture 1_b");
    hold off;
    
    % subplot(2,2,2);
    % hold on;
    % yyaxis left;
    % plot(t_res,y_res{4},'b-','LineWidth',1.5);
    % ylabel("Concentration (g/L)");
    % axes = fig1_b.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % plot(t_res,y_res{6},'k-','LineWidth',1.5);
    % legend(["Glucose","Product"]);
    % xlabel("Time (hr)"); ylabel("Concentration (g/L)");
    % title("Repeating Batch Culture");
    % hold off;
    % 
    % subplot(2,2,3);
    % hold on;
    % plot(t_res,y_res{5},'b-','LineWidth',1.5);
    % xlabel("Time (hr)"); ylabel("Dissolved O_2 Concentration (g/L)");
    % ylim([0,0.15])
    % hold off;
    % 
    % subplot(2,2,4);
    % hold on;
    % plot(0:length(tot_prod_mass_cum)-1,tot_prod_mass_cum,'b-','LineWidth',1.5);
    % xlabel("Batch Number"); ylabel("Total Product Mass (g)");
    % hold off;

    axes = fig1_b.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err
    
        end
    end

    disp(table(tot_prod_mass_cum(end).*1000,'VariableNames',{'Batch'}))

    % 3 different strains - O2 oscillation
    % Environmental Conditions
    dV = 0.06; % L/s
    V = 1.5;% L
    C_G_in = 20; % g/L
    S_O2 = 0.007;
    k_L_a = 10;
    % k_L_a = 2;
    amp = 0.075;
    % amp = 0;
    freq = 50;
    y_shift = 0.3;
    env_cond = [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift];

    fig1_e = figure;
    fig1_e.Color = [1,1,1];

    ODE_4 = ode;
    ODE_4.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
    ODE_4.InitialTime = 0;
    ODE_4.InitialValue = [1.25,0.25,0.05,20,0.007,0,0];
    ODE_4.Parameters = p;
    ODE_4.Solver = 'ode15s';
    
    sol = solve(ODE_4,t_span(1),t_span(2));
    t_res_1_e = sol.Time;
    y_res_1_e = sol.Solution;

    hold on;
    % yyaxis left;
    plot(t_res_1_e,y_res_1_e(1,:),'r-','LineWidth',1.5);
    plot(t_res_1_e,y_res_1_e(2,:),'k-','LineWidth',1.5);
    plot(t_res_1_e,y_res_1_e(3,:),'b-','LineWidth',1.5);
    plot(t_res_1_e,sum(y_res_1_e(1:3,:),1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
    axes = fig1_e.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    % yyaxis right;
    % plot(t_res_1_e,y_res_1_e(6,:).*1000,'b--','LineWidth',1.5);
    % ylim([0,0.15.*1000]);
    % legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass"],'Location','best');
    xlabel("Time (hr)"); % ylabel("\beta-Carotene Concentration (mg/L)")
    title("O2 oscillation 1_e");
    hold off;
    
    axes = fig1_e.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end

    fig1_f = figure;
    fig1_f.Color = [1,1,1];
    hold on;
    plot(t_res_1_e,y_res_1_e(5,:),'b-','LineWidth',1.5);
    xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
    title("O2 oscillation 1_f")
    ylim([0,Inf]);
    hold off;

    axes = fig1_f.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end

    % 3 different strains - O2 limitation
    % Environmental Conditions
    dV = 0.06; % L/s
    V = 1.5;% L
    C_G_in = 20; % g/L
    S_O2 = 0.007;
    % k_L_a = 10;
    k_L_a = 2;
    % amp = 0.075;
    amp = 0;
    freq = 50;
    y_shift = 0.3;
    env_cond = [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift];

    fig1_g = figure;
    fig1_g.Color = [1,1,1];

    ODE_4 = ode;
    ODE_4.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
    ODE_4.InitialTime = 0;
    ODE_4.InitialValue = [1.25,0.25,0.05,20,0.007,0,0];
    ODE_4.Parameters = p;
    ODE_4.Solver = 'ode15s';
    
    sol = solve(ODE_4,t_span(1),t_span(2));
    t_res_1_g = sol.Time;
    y_res_1_g = sol.Solution;
    
    hold on;
    % yyaxis left;
    plot(t_res_1_g,y_res_1_g(1,:),'r-','LineWidth',1.5);
    plot(t_res_1_g,y_res_1_g(2,:),'k-','LineWidth',1.5);
    plot(t_res_1_g,y_res_1_g(3,:),'b-','LineWidth',1.5);
    plot(t_res_1_g,sum(y_res_1_g(1:3,:),1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
    axes = fig1_g.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    % yyaxis right;
    % plot(t_res_1_g,y_res_1_g(6,:).*1000,'b--','LineWidth',1.5);
    % ylim([0,0.15.*1000]);
    % legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass"],'Location','best');
    xlabel("Time (hr)"); % ylabel("\beta-Carotene Concentration (mg/L)")
    title("Low kLa 1_g");
    hold off;

    axes = fig1_g.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    
    fig1_h = figure;
    fig1_h.Color = [1,1,1];
    hold on;
    plot(t_res_1_g,y_res_1_g(5,:),'b-','LineWidth',1.5);
    xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
    title("Low kLa 1_h")
    ylim([0,Inf]);
    hold off;

    axes = fig1_h.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end

    disp(table(y_res_1_e(7,end).*1000,y_res_1_g(7,end).*1000,'VariableNames',{'O2 oscillation','low kLa'}))

    % combined production plot 1
    fig_last = figure;
    fig_last.Color = [1,1,1];
    hold on;
    % 1_a production
    plot(t_res_1_a,y_res_1_a(6,:).*1000,'r-','LineWidth',1.5);
    % 1_b production
    plot(t_res_1_b,y_res_1_b{6}.*1000,'k--','LineWidth',1);
    % 1_d production
    plot(t_res_1_d,y_res_1_d(6,:).*1000,'g-','LineWidth',1.5);
    % 1_e production
    plot(t_res_1_e,y_res_1_e(6,:).*1000,'b-','LineWidth',1.5);
    % 1_g production
    plot(t_res_1_g,y_res_1_g(6,:).*1000,'k-','LineWidth',1.5);
    ylim([0,0.15.*1000]);
    legend(["Fig. 1A","Fig. 1B","Fig. 1C","Fig. 1D","Fig. 1E"]);
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)");
    hold off;
    
    axes = fig_last.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end

    % manuscript plots
    t_span = [0,1000];

    C_G_in_new = 20; % g/L

    dilution_rate_factor = 2;
    dV_range = [0.06./dilution_rate_factor,0.06,0.06.*dilution_rate_factor];

    plot_style_var = {'-','--','o-'};
    plot_color_var = {'r','k','b'};

    % a - higher dilution rate => more 
    % fig1_z = figure;
    % fig1_z.Color = [1,1,1];
    % 
    % for k=1:1:length(dV_range)
    %     % Environmental Conditions
    %     dV = dV_range(k); % L/s
    %     V = 1.5;% L
    %     S_O2 = 0.007;
    %     k_L_a = 10;
    %     amp = 0;
    %     freq = 50;
    %     y_shift = 0.3;
    %     env_cond = [dV,V,C_G_in_new,S_O2,k_L_a,amp,freq,y_shift];
    % 
    %     alpha_HI_prod = 0.001.*alpha_mult;
    %     alpha_MED_prod = 0.00075.*alpha_mult;
    %     alpha_LO_prod = 0.0005.*alpha_mult;
    % 
    %     beta_HI_prod = 0.001.*beta_mult.*0;
    %     beta_MED_prod = 0.00075.*beta_mult.*0;
    %     beta_LO_prod = 0.0005.*beta_mult.*0;
    % 
    %     p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
    %         mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
    %         alpha_HI_prod,alpha_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,alpha_MED_prod, ...
    %         beta_HI_prod,beta_MED_prod,beta_LO_prod,Y_prod_G];
    % 
    %     ODE_7 = ode;
    %     ODE_7.ODEFcn = @(t,y,p) ODESys_cont_3_alpha_beta_var(t,y,p,env_cond);
    %     ODE_7.InitialTime = 0;
    %     ODE_7.InitialValue = [1.25,0.25,0.05,20,0.007,0,0,0,0];
    %     ODE_7.Parameters = p;
    %     ODE_7.Solver = 'ode15s';
    % 
    %     sol = solve(ODE_7,t_span(1),t_span(2));
    %     t_res_1_z = sol.Time;
    %     y_res_1_z = sol.Solution;
    % 
    %     yyaxis left;
    %     hold on;
    %     plot(t_res_1_z,y_res_1_z(6,:).*1000,[plot_color_var{k},'-'],'LineWidth',1.5);
    %     ylabel("\beta-Carotene Concentration (mg/L)");
    %     hold off;
    %     yyaxis right;
    %     hold on;
    %     plot(t_res_1_z,y_res_1_z(7,:).*1000,[plot_color_var{k},'--'],'LineWidth',1.5);
    %     ylabel("Total \beta-Carotene Mass (mg)");
    %     hold off;
    % end
    % 
    % legend(["low","med","high","low total","med total","high total"]);
    % title("z");
    % yyaxis left;
    % ylim([0,150]);
    % axes = fig1_z.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % ylim([0,3000]);
    % axes = fig1_z.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % if LO_multiplier == 0.95
    %     if C_G_in_new == 20
    %         saveas(fig1_z,'fig_a.png');
    %         saveas(fig1_z,'fig_a.fig');
    %     end
    % end
    % 
    % % b - low vs high glucose feed concentration
    % fig1_y = figure;
    % fig1_y.Color = [1,1,1];
    % 
    % glucose_feed = [20,50]; % g/L
    % for k=1:1:length(glucose_feed)
    %     % Environmental Conditions
    %     dV = 0.06; % L/s
    %     V = 1.5;% L
    %     C_G_in_3 = glucose_feed(k);
    %     S_O2 = 0.007;
    %     k_L_a = 10;
    %     amp = 0;
    %     freq = 50;
    %     y_shift = 0.3;
    %     env_cond = [dV,V,C_G_in_new,S_O2,k_L_a,amp,freq,y_shift];
    % 
    %     alpha_HI_prod = 0.001.*alpha_mult;
    %     alpha_MED_prod = 0.00075.*alpha_mult;
    %     alpha_LO_prod = 0.0005.*alpha_mult;
    % 
    %     beta_HI_prod = 0.001.*beta_mult;
    %     beta_MED_prod = 0.00075.*beta_mult;
    %     beta_LO_prod = 0.0005.*beta_mult;
    % 
    %     p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
    %         mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
    %         alpha_HI_prod,alpha_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,alpha_MED_prod, ...
    %         beta_HI_prod,beta_MED_prod,beta_LO_prod,Y_prod_G];
    % 
    %     ODE_7 = ode;
    %     ODE_7.ODEFcn = @(t,y,p) ODESys_cont_3_alpha_beta_var(t,y,p,env_cond);
    %     ODE_7.InitialTime = 0;
    %     ODE_7.InitialValue = [1.25,0.25,0.05,20,0.007,0,0,0,0];
    %     ODE_7.Parameters = p;
    %     ODE_7.Solver = 'ode15s';
    % 
    %     sol = solve(ODE_7,t_span(1),t_span(2));
    %     t_res_1_y = sol.Time;
    %     y_res_1_y = sol.Solution;
    % 
    %     yyaxis left;
    %     hold on;
    %     plot(t_res_1_y,y_res_1_y(6,:).*1000,[plot_color_var{k},'-'],'LineWidth',1.5);
    %     ylabel("\beta-Carotene Concentration (mg/L)");
    %     hold off;
    %     yyaxis right;
    %     hold on;
    %     plot(t_res_1_y,y_res_1_y(7,:).*1000,[plot_color_var{k},'--'],'LineWidth',1.5);
    %     ylabel("Total \beta-Carotene Mass (mg)");
    %     hold off;
    % end
    % 
    % title("y");
    % yyaxis left;
    % ylim([0,150]);
    % axes = fig1_y.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % yyaxis right;
    % ylim([0,3000]);
    % axes = fig1_y.Children;
    % for k=1:1:length(axes)
    %     axes(k).FontName = 'Arial';
    %     axes(k).FontSize = 15;
    %     axes(k).FontWeight = 'bold';
    %     axes(k).Box = true;
    %     try
    %         axes(k).YColor = [0,0,0];
    %     catch err
    % 
    %     end
    % end
    % if LO_multiplier == 0.8
    %     if C_G_in_new == 20
    %         saveas(fig1_y,'fig_b.png');
    %         saveas(fig1_y,'fig_b.fig');
    %     end
    % end

    % 3 different strains - alpha only - product
    fig1_i = figure;
    fig1_i.Color = [1,1,1];

    for k=1:1:length(dV_range)
        % Environmental Conditions
        dV = dV_range(k); % L/s
        V = 1.5;% L
        S_O2 = 0.007;
        k_L_a = 10;
        amp = 0;
        freq = 50;
        y_shift = 0.3;
        env_cond = [dV,V,C_G_in_new,S_O2,k_L_a,amp,freq,y_shift];

        alpha_HI_prod = 0.001.*alpha_mult;
        alpha_MED_prod = 0.00075.*alpha_mult;
        alpha_LO_prod = 0.0005.*alpha_mult;
    
        beta_HI_prod = 0.001.*beta_mult.*0;
        beta_MED_prod = 0.00075.*beta_mult.*0;
        beta_LO_prod = 0.0005.*beta_mult.*0;
    
        p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
            mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
            alpha_HI_prod,alpha_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,alpha_MED_prod, ...
            beta_HI_prod,beta_MED_prod,beta_LO_prod,Y_prod_G];

        ODE_7 = ode;
        ODE_7.ODEFcn = @(t,y,p) ODESys_cont_3_alpha_beta_var(t,y,p,env_cond);
        ODE_7.InitialTime = 0;
        ODE_7.InitialValue = [1.25,0.25,0.05,20,0.007,0,0,0,0];
        ODE_7.Parameters = p;
        ODE_7.Solver = 'ode15s';
        
        sol = solve(ODE_7,t_span(1),t_span(2));
        t_res_1_i = sol.Time;
        y_res_1_i = sol.Solution;

        yyaxis left;
        hold on;
        plot(t_res_1_i,y_res_1_i(6,:).*1000,[plot_color_var{k},'-'],'LineWidth',1.5);
        ylabel("\beta-Carotene Concentration (mg/L)");
        hold off;
        yyaxis right;
        hold on;
        plot(t_res_1_i,y_res_1_i(7,:).*1000,[plot_color_var{k},'--'],'LineWidth',1.5);
        ylabel("Total \beta-Carotene Mass (mg)");
        hold off;
    end

    % legend(["low","med","high","low total","med total","high total"]);
    yyaxis left;
    ylim([0,500]);
    axes = fig1_i.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    yyaxis right;
    ylim([0,10000]);
    axes = fig1_i.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    if LO_multiplier == 0.8
        if C_G_in_new == 20
            saveas(fig1_i,'alpha_lowComp_lowGlu.png')
            saveas(fig1_i,'alpha_lowComp_lowGlu.fig')
        else
            saveas(fig1_i,'alpha_lowComp_hiGlu.png')
            saveas(fig1_i,'alpha_lowComp_hiGlu.fig')
        end
    else
        if C_G_in_new == 20
            saveas(fig1_i,'alpha_hiComp_lowGlu.png')
            saveas(fig1_i,'alpha_hiComp_lowGlu.fig')
        else
            saveas(fig1_i,'alpha_hiComp_hiGlu.png')
            saveas(fig1_i,'alpha_hiComp_hiGlu.fig')
        end
    end

    % 3 different strains - alpha only - biomass
    fig1_k = figure;
    fig1_k.Color = [1,1,1];

    for k=1:1:length(dV_range)
        % Environmental Conditions
        dV = dV_range(k); % L/s
        V = 1.5; % L
        S_O2 = 0.007;
        k_L_a = 10;
        amp = 0;
        freq = 50;
        y_shift = 0.3;
        env_cond = [dV,V,C_G_in_new,S_O2,k_L_a,amp,freq,y_shift];

        alpha_HI_prod = 0.001.*alpha_mult;
        alpha_MED_prod = 0.00075.*alpha_mult;
        alpha_LO_prod = 0.0005.*alpha_mult;
    
        beta_HI_prod = 0.001.*beta_mult;
        beta_MED_prod = 0.00075.*beta_mult;
        beta_LO_prod = 0.0005.*beta_mult;
    
        p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
            mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
            alpha_HI_prod,alpha_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,alpha_MED_prod, ...
            beta_HI_prod,beta_MED_prod,beta_LO_prod,Y_prod_G];
    
        ODE_8 = ode;
        ODE_8.ODEFcn = @(t,y,p) ODESys_cont_3_alpha_beta_var(t,y,p,env_cond);
        ODE_8.InitialTime = 0;
        ODE_8.InitialValue = [1.25,0.25,0.05,20,0.007,0,0,0,0];
        ODE_8.Parameters = p;
        ODE_8.Solver = 'ode15s';
        
        sol = solve(ODE_8,t_span(1),t_span(2));
        t_res_1_k = sol.Time;
        y_res_1_k = sol.Solution;
        
        hold on;
        yyaxis left;
        plot(t_res_1_k,y_res_1_k(1,:),['r',plot_style_var{k}],'LineWidth',1.5);
        plot(t_res_1_k,y_res_1_k(2,:),['k',plot_style_var{k}],'LineWidth',1.5);
        plot(t_res_1_k,y_res_1_k(3,:),['b',plot_style_var{k}],'LineWidth',1.5);
        ylabel("Biomass (g/L)");
        yyaxis right;
        plot(t_res_1_k,y_res_1_k(4,:),['g',plot_style_var{k}],'LineWidth',1.5);      
        ylabel("Glucose (g/L)");
        hold off;
    end

    yyaxis left;
    ylim([0,20]);
    axes = fig1_k.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    yyaxis right;
    ylim([0,40]);
    axes = fig1_k.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    hold off;

    if LO_multiplier == 0.8
        if C_G_in_new == 20
            saveas(fig1_k,'biomass_lowComp_lowGlu.png')
            saveas(fig1_k,'biomass_lowComp_lowGlu.fig')
        else
            saveas(fig1_k,'biomass_lowComp_hiGlu.png')
            saveas(fig1_k,'biomass_lowComp_hiGlu.fig')
        end
    else
        if C_G_in_new == 20
            saveas(fig1_k,'biomass_hiComp_lowGlu.png')
            saveas(fig1_k,'biomass_hiComp_lowGlu.fig')
        else
            saveas(fig1_k,'biomass_hiComp_hiGlu.png')
            saveas(fig1_k,'biomass_hiComp_hiGlu.fig')
        end
    end

    % 3 different strains - beta only - product
    fig1_j = figure;
    fig1_j.Color = [1,1,1];

    for k=1:1:length(dV_range)
        % Environmental Conditions
        dV = dV_range(k); % L/s
        V = 1.5; % L
        S_O2 = 0.007;
        k_L_a = 10;
        amp = 0;
        freq = 50;
        y_shift = 0.3;
        env_cond = [dV,V,C_G_in_new,S_O2,k_L_a,amp,freq,y_shift];

        alpha_HI_prod = 0.001.*alpha_mult.*0;
        alpha_MED_prod = 0.00075.*alpha_mult.*0;
        alpha_LO_prod = 0.0005.*alpha_mult.*0;
    
        beta_HI_prod = 0.001.*beta_mult;
        beta_MED_prod = 0.00075.*beta_mult;
        beta_LO_prod = 0.0005.*beta_mult;
    
        p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
            mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
            alpha_HI_prod,alpha_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,alpha_MED_prod, ...
            beta_HI_prod,beta_MED_prod,beta_LO_prod,Y_prod_G];
    
        ODE_8 = ode;
        ODE_8.ODEFcn = @(t,y,p) ODESys_cont_3_alpha_beta_var(t,y,p,env_cond);
        ODE_8.InitialTime = 0;
        ODE_8.InitialValue = [1.25,0.25,0.05,20,0.007,0,0,0,0];
        ODE_8.Parameters = p;
        ODE_8.Solver = 'ode15s';
        
        sol = solve(ODE_8,t_span(1),t_span(2));
        t_res_1_j = sol.Time;
        y_res_1_j = sol.Solution;
        
        yyaxis left;
        hold on;
        plot(t_res_1_j,y_res_1_j(6,:).*1000,[plot_color_var{k},'-'],'LineWidth',1.5);
        ylabel("\beta-Carotene Concentration (mg/L)");
        hold off;
        yyaxis right;
        hold on;
        plot(t_res_1_j,y_res_1_j(7,:).*1000,[plot_color_var{k},'--'],'LineWidth',1.5);
        ylabel("Total \beta-Carotene Mass (mg)");
        hold off;
    end

    yyaxis left;
    ylim([0,500]);
    axes = fig1_j.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    yyaxis right;
    ylim([0,10000]);
    axes = fig1_j.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    hold off;

    if LO_multiplier == 0.8
        if C_G_in_new == 20
            saveas(fig1_j,'beta_lowComp_lowGlu.png')
            saveas(fig1_j,'beta_lowComp_lowGlu.fig')
        else
            saveas(fig1_j,'beta_lowComp_hiGlu.png')
            saveas(fig1_j,'beta_lowComp_hiGlu.fig')
        end
    else
        if C_G_in_new == 20
            saveas(fig1_j,'beta_hiComp_lowGlu.png')
            saveas(fig1_j,'beta_hiComp_lowGlu.fig')
        else
            saveas(fig1_j,'beta_hiComp_hiGlu.png')
            saveas(fig1_j,'beta_hiComp_hiGlu.fig')
        end
    end

    % 3 different strains - balanced alpha beta check
    fig1_m = figure;
    fig1_m.Color = [1,1,1];

    % Environmental Conditions
    dV = 0.06; % L/s
    V = 1.5;% L
    S_O2 = 0.007;
    k_L_a = 10;
    amp = 0;
    freq = 50;
    y_shift = 0.3;
    env_cond = [dV,V,C_G_in_new,S_O2,k_L_a,amp,freq,y_shift];

    alpha_HI_prod = 0.001.*alpha_mult;
    alpha_MED_prod = 0.00075.*alpha_mult;
    alpha_LO_prod = 0.0005.*alpha_mult;

    beta_HI_prod = 0.001.*beta_mult;
    beta_MED_prod = 0.00075.*beta_mult;
    beta_LO_prod = 0.0005.*beta_mult;

    p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
        mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
        alpha_HI_prod,alpha_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,alpha_MED_prod, ...
        beta_HI_prod,beta_MED_prod,beta_LO_prod,Y_prod_G];

    ODE_9 = ode;
    ODE_9.ODEFcn = @(t,y,p) ODESys_cont_3_alpha_beta_var(t,y,p,env_cond);
    ODE_9.InitialTime = 0;
    ODE_9.InitialValue = [1.25,0.25,0.05,20,0.007,0,0,0,0];
    ODE_9.Parameters = p;
    ODE_9.Solver = 'ode15s';

    sol = solve(ODE_9,t_span(1),t_span(2));
    t_res_1_l = sol.Time;
    y_res_1_l = sol.Solution;

    hold on;
    plot(t_res_1_l,y_res_1_l(8,:).*1000,'LineWidth',1.5);
    plot(t_res_1_l,y_res_1_l(9,:).*1000,'LineWidth',1.5);
    plot(t_res_1_l,y_res_1_l(6,:).*1000,'LineWidth',1.5);
    hold off;

    ylabel("\beta-Carotene Concentration (mg/L)");
    legend(["alpha","beta","tot"]);
    title("balance");
    ylim([0,max(y_res_1_l([6,8,9],:),[],"All").*1000]);
    axes = fig1_m.Children;
    for k=1:1:length(axes)
        axes(k).FontName = 'Arial';
        axes(k).FontSize = 15;
        axes(k).FontWeight = 'bold';
        axes(k).Box = true;
        try
            axes(k).YColor = [0,0,0];
        catch err

        end
    end
    hold off;
    
else
    disp('Error: Please set MED_multiplier and LO_multipler within the range [0.6,1].')
end

%% Functions
function v = EvtFcn(t,y,p)
    v(1) = y(4) - 0.1;
end

function dydt = simple_mu_exper(t,y,p)
    dV = 0.06;
    V = 1.5;
    C_in = 20;
    dydt(1) = y(1).*p(1).*(y(3)./(y(3)+p(3))) - dV./V.*y(1);
    dydt(2) = y(2).*p(2).*(y(3)./(y(3)+p(4))) - dV./V.*y(2);
    dydt(3) = -p(5).*y(1).*p(1).*(y(3)./(y(3)+p(3))) - p(6).*y(2).*p(2).*(y(3)./(y(3)+p(4))) + dV./V.*C_in - dV./V.*y(3);

    dydt = dydt';
end

function dydt = ODESys_cont_2(t,y,p,env_cond)
    % 1. HI
    % 2. LO
    % 3. NO
    % 4. Glucose
    % 5. O2
    % 6. prod
    % 7. prod mass
    
    env_cond_cell = num2cell(env_cond);
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift] = deal(env_cond_cell{:});
    % if t < 24, dV = 0; end

    dydt(1) = p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(4).*y(1) - ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - dV./V.*y(1);
    % dydt(2) = p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - p(9).*y(2) - ...
    %     p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - dV./V.*y(2);
    
    dydt(2) = 0;
    
    dydt(3) = p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(14).*y(3) + ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - dV./V.*y(3);
    
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) ...
        - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) ...
        - p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) ...
        - p(34).*y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        - p(34).*y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)))) ...
        + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(pi./24)*(cos(t./24.*pi));
    
    dydt(6) = y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        + y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        - dV./V.*y(6);
    
    dydt(7) = dV.*y(6);

    dydt = dydt';
end

function dydt = ODESys_batch_3(t,y,p,env_cond)
    % 1. HI
    % 2. LO
    % 3. NO
    % 4. Glucose
    % 5. O2
    % 6. prod
    % 7. prod mass
    
    env_cond_cell = num2cell(env_cond);
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift] = deal(env_cond_cell{:});
    % if t < 24, dV = 0; end

    dydt(1) = p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(4).*y(1) - ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)));
    
    dydt(2) = p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - p(9).*y(2) - ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)));
    
    dydt(3) = p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(14).*y(3) + ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)));
    
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) ...
        - p(34).*y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        - p(34).*y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)))) ...
        - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5)));
    % dydt(4) = 0;
    
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(cos(freq.*(2./pi).*t));
    
    dydt(6) = y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        + y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))));
    % dydt(7) = dV.*y(6);

    dydt = dydt';
end

function dydt = ODESys_cont_3_no_mut(t,y,p,env_cond)
    % 1. HI
    % 2. LO
    % 3. NO
    % 4. Glucose
    % 5. O2
    % 6. prod
    % 7. prod mass
    
    env_cond_cell = num2cell(env_cond);
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift] = deal(env_cond_cell{:});
    % if t < 24, dV = 0; end

    dydt(1) = p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(4).*y(1) - ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - dV./V.*y(1);
    
    dydt(2) = p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - p(9).*y(2) - ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - dV./V.*y(2);
    
    dydt(3) = p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(14).*y(3) - dV./V.*y(3);
        % p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + ...
        % p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))
    
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) ...
        - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) ...
        - p(34).*y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        - p(34).*y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)))) ...
        - p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) ...
        + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(cos(freq.*(2./pi).*t));
    
    dydt(6) = y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        + y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)))) ...
        - dV./V.*y(6);
    
    dydt(7) = dV.*y(6);

    dydt = dydt';
end

function dydt = ODESys_cont_3(t,y,p,env_cond)
    % 1. HI
    % 2. LO
    % 3. NO
    % 4. Glucose
    % 5. O2
    % 6. prod
    % 7. prod mass

    env_cond_cell = num2cell(env_cond);
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift] = deal(env_cond_cell{:});
    % if t < 24, dV = 0; end

    dydt(1) = p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(4).*y(1) - ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - dV./V.*y(1);

    dydt(2) = p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - p(9).*y(2) - ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - dV./V.*y(2);

    dydt(3) = p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(14).*y(3) + ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - dV./V.*y(3);
    
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) ...
        - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) ...
        - p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) ...
        - p(34).*y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        - p(34).*y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)))) ...
        + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(cos(freq.*(2./pi).*t));
    
    dydt(6) = y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        + y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)))) - ...
        dV./V.*y(6);
    
    dydt(7) = dV.*y(6);

    dydt = dydt';
end

function dydt = ODESys_cont_3_alpha_beta_var(t,y,p,env_cond)
    % 1. HI
    % 2. LO
    % 3. NO
    % 4. Glucose
    % 5. O2
    % 6. prod
    % 7. prod mass
    % 8. prod_growth_associated
    % 9. prod_non_growth_associated

    max_biomass = 20; % g/L

    env_cond_cell = num2cell(env_cond);
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift] = deal(env_cond_cell{:});
    % if t < 24, dV = 0; end

    dydt(1) = p(1).*y(1).*(y(4)./(p(2)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) - p(4).*y(1) - ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) - dV./V.*y(1);
    
    dydt(2) = p(6).*y(2).*(y(4)./(p(7)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) - p(9).*y(2) - ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) - dV./V.*y(2);
    
    dydt(3) = p(11).*y(3).*(y(4)./(p(12)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) - p(14).*y(3) + ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) + ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) - dV./V.*y(3);
    
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) ...
        - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) ...
        - p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) ...
        - p(34).*y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) + p(31).*(y(4)./(p(2)+y(4)))) ...
        - p(34).*y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) + p(33).*(y(4)./(p(7)+y(4)))) ...
        + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    
    % dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
    %     p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(cos(freq.*(2./pi).*t));
    dydt(5) = 0;
    
    dydt(6) = y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) + p(31).*(y(4)./(p(2)+y(4)))) ...
        + y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) + p(33).*(y(4)./(p(7)+y(4)))) ...
        - dV./V.*y(6);
    
    dydt(7) = dV.*y(6);
    
    dydt(8) = y(1).*p(21).*p(1).*(y(4)./(p(2)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) ...
        + y(2).*p(22).*p(6).*(y(4)./(p(7)+y(4))).*(1 - (y(1)+y(2)+y(3))./max_biomass) ...
        - dV./V.*y(8);
    
    dydt(9) = y(1).*p(31).*(y(4)./(p(2)+y(4))) + y(2).*p(33).*(y(4)./(p(7)+y(4))) - dV./V.*y(9);

    dydt = dydt';
end

function dydt = ODESys_cont_4(t,y,p,env_cond)
    % 1. HI
    % 2. LO
    % 3. NO
    % 4. Glucose
    % 5. O2
    % 6. prod
    % 7. prod mass
    % 8. MED
        
    env_cond_cell = num2cell(env_cond);
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq,y_shift] = deal(env_cond_cell{:});
    % if t < 24, dV = 0; end

    dydt(1) = p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(4).*y(1) - ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - dV./V.*y(1);
    
    dydt(2) = p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - p(9).*y(2) - ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - dV./V.*y(2);
    
    dydt(3) = p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(14).*y(3) + ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + ...
        p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + ...
        p(27).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - dV./V.*y(3);
    
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) ...
        - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) ...
        - p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) ...
        - p(28).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) ...
        - p(34).*y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        - p(34).*y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)))) ...
        - p(34).*y(8).*(p(30).*p(23).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) + p(32).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5)))) ...
        + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(29).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(pi./24)*(cos(t./24.*pi));
    
    dydt(6) = y(1).*(p(21).*p(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(31).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)))) ...
        + y(2).*(p(22).*p(6).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(33).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)))) ...
        + y(8).*(p(30).*p(23).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) + p(32).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5)))) ...
        - dV./V.*y(6);
    
    dydt(7) = dV.*y(6);
    
    dydt(8) = p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - p(26).*y(8) ...
        - p(27).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - dV./V.*y(8);
    % dydt(8) = 0;

    dydt = dydt';
end