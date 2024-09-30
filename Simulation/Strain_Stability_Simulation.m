clear;
clc;
close all;

%% User Input
% Modify the maximum specific growth rates of the intermediate strains
% (medium-producing [MED] and low-producing [LO]) values relative to that
% of the wild type (WT).

% MED maximum specific growth rate multiplier (valid values: [0.6,1]):
MED_multiplier = 0.65;

% LO maximum specific growth rate multipler (valid values: [0.6,1]):
LO_multiplier = 0.7;

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
    
    Y_HI_prod = 0.1;
    Y_MED_prod = 0.075;
    Y_LO_prod = 0.05;
    
    mu_max_4 = 0.3466.*MED_multiplier;
    K_MED_G = 1;
    K_MED_O2 = 0.002./MED_multiplier;
    k_d_MED = 0.005;
    MR_MED = 0.00005;
    
    p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
            mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
            Y_HI_prod,Y_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,Y_MED_prod];
    
    param_names = "mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO,mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2,Y_HI_prod,Y_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,Y_MED_prod";
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
    
    % base dilution rate (0.05 1/hr)
    ODE_1 = ode;
    ODE_1.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
    ODE_1.InitialTime = 0;
    ODE_1.InitialValue = [1.25,0.25,0.05,20,0.007,0,0];
    ODE_1.Parameters = p;
    ODE_1.Solver = 'ode15s';
    
    sol = solve(ODE_1,t_span(1),t_span(2));
    t_res = sol.Time;
    y_res = sol.Solution;
    
    % subplot(1,2,1);
    hold on;
    yyaxis left;
    plot(t_res,y_res(1,:),'b-','LineWidth',1.5);
    plot(t_res,y_res(2,:),'k-','LineWidth',1.5);
    plot(t_res,y_res(3,:),'r-','LineWidth',1.5);
    plot(t_res,sum(y_res(1:3,:),1),'k--','LineWidth',1.5);
    ylabel("Biomass Concentration (g/L)");
    ylim([0,10]);
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
    yyaxis right;
    plot(t_res,y_res(6,:).*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000]);
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)");
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
    
    % higher dilution rate (0.1 1/hr)
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
    t_res_2 = sol.Time;
    y_res_2 = sol.Solution;
    
    % subplot(1,2,2);
    hold on;
    plot(t_res_2,y_res_2(1,:),'b-','LineWidth',1.5);
    plot(t_res_2,y_res_2(2,:),'k-','LineWidth',1.5);
    plot(t_res_2,y_res_2(3,:),'r-','LineWidth',1.5);
    plot(t_res_2,sum(y_res_2(1:3,:),1),'k--','LineWidth',1.5);
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
    yyaxis right;
    plot(t_res_2,y_res_2(6,:).*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000]);
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)");
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

    disp(table(y_res(7,end).*1000,y_res_2(7,end).*1000,'VariableNames',{'0.05 1/hr','0.1 1/hr'}))
    
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
    t_res = sol.Time;
    y_res = sol.Solution;
    
    hold on;
    yyaxis left;
    plot(t_res,y_res(1,:),'b-','LineWidth',1.5);
    plot(t_res,y_res(3,:),'r-','LineWidth',1.5);
    plot(t_res,sum(y_res(1:3,:),1),'k--','LineWidth',1.5);
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
    plot(t_res,y_res(6,:).*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000]);
    legend(["High-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
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
    t_res_2 = sol.Time;
    y_res_2 = sol.Solution;
    
    subplot(1,2,1);
    hold on;
    yyaxis left;
    plot(t_res_2,y_res_2(1,:),'b-','LineWidth',1.5);
    plot(t_res_2,y_res_2(2,:),'k-','LineWidth',1.5);
    plot(t_res_2,y_res_2(3,:),'r-','LineWidth',1.5);
    plot(t_res_2,sum(y_res_2(1:3,:),1),'k--','LineWidth',1.5);
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
    plot(t_res_2,y_res_2(6,:).*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000]);
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
    t_res_3 = sol.Time;
    y_res_3 = sol.Solution;
    
    subplot(1,2,2);
    hold on;
    yyaxis left;
    plot(t_res_3,y_res_3(1,:),'b-','LineWidth',1.5);
    plot(t_res_3,y_res_3(8,:),'g-','LineWidth',1.5);
    plot(t_res_3,y_res_3(2,:),'k-','LineWidth',1.5);
    plot(t_res_3,y_res_3(3,:),'r-','LineWidth',1.5);
    plot(t_res_3,sum([y_res_3(1:3,:);y_res_3(8,:)],1),'k--','LineWidth',1.5);
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
    plot(t_res_3,y_res_3(6,:).*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000]);
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
    yyaxis left;
    plot(t_res_4,y_res_4(1,:),'b-','LineWidth',1.5);
    plot(t_res_4,y_res_4(2,:),'k-','LineWidth',1.5);
    plot(t_res_4,y_res_4(3,:),'r-','LineWidth',1.5);
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
    plot(t_res_4,y_res_4(6,:).*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000]);
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
    
    disp(table(y_res(7,end).*1000,y_res_2(7,end).*1000,y_res_3(7,end).*1000,y_res_4(7,end).*1000,'VariableNames',{'2','3','4','3 - no mut'}))
    
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
    t_res = [];
    y_res = cell(6,1);
    tot_prod_mass = 0;
    V = 0.05; % L
    batch_time_length = 24; % hr
    t_span_1 = [0,batch_time_length];
    ODE_7 = ode;
    ODE_7.ODEFcn = @(t,y,p) ODESys_batch_3(t,y,p,env_cond);
    ODE_7.InitialTime = 0;
    ODE_7.InitialValue = [0.125,0.025,0,20,0.007,0];
    ODE_7.Parameters = p;
    ODE_7.Solver = 'ode15s';
    % evt = odeEvent(EventFcn=@(t,y,p) EvtFcn(t,y,p),Response="stop");
    % ODE_7.EventDefinition = evt;
    while (batch_ct < 9)
        sol = solve(ODE_7,t_span_1(1),t_span_1(2));

        t_res = [t_res;sol.Time']; %#ok<*AGROW>
        for l=1:1:size(sol.Solution,1)
            y_res{l} = [y_res{l},sol.Solution(l,:)]; %#ok<*SAGROW>
        end

        t_span_1(1) = sol.Time(end);
        ODE_7.InitialTime = t_span_1(1);
        t_span_1(2) = t_span_1(1)+batch_time_length;

        ODE_7.InitialValue = [y_res{1}(end).*0.01,y_res{2}(end).*0.01,y_res{3}(end).*0.01,20,0.007,0];
        % ODE_7.InitialValue = [0.125,0.025,0.005,20,0.007,0];

        tot_prod_mass(end+1) = y_res{6}(end).*V;

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
    yyaxis left;
    plot(t_res,y_res{1},'b-','LineWidth',1.5);
    plot(t_res,y_res{2},'k-','LineWidth',1.5);
    plot(t_res,y_res{3},'r-','LineWidth',1.5);
    plot(t_res,sum([y_res{1};y_res{2};y_res{3}],1),'k--','LineWidth',1.5);
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
    yyaxis right;
    plot(t_res,y_res{6}.*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000])
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"]);
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)");
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
    % ylim([0,0.1])
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
    yyaxis left;
    plot(t_res_1_e,y_res_1_e(1,:),'b-','LineWidth',1.5);
    plot(t_res_1_e,y_res_1_e(2,:),'k-','LineWidth',1.5);
    plot(t_res_1_e,y_res_1_e(3,:),'r-','LineWidth',1.5);
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
    yyaxis right;
    plot(t_res_1_e,y_res_1_e(6,:).*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000]);
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)")
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

    % 3 different strains - O2 oscillation
    % Environmental Conditions
    dV = 0.06; % L/s
    V = 1.5;% L
    C_G_in = 20; % g/L
    S_O2 = 0.007;
    % k_L_a = 10;
    k_L_a = 2;
    % pro0.075;
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
    yyaxis left;
    plot(t_res_1_g,y_res_1_g(1,:),'b-','LineWidth',1.5);
    plot(t_res_1_g,y_res_1_g(2,:),'k-','LineWidth',1.5);
    plot(t_res_1_g,y_res_1_g(3,:),'r-','LineWidth',1.5);
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
    yyaxis right;
    plot(t_res_1_g,y_res_1_g(6,:).*1000,'b--','LineWidth',1.5);
    ylim([0,0.1.*1000]);
    legend(["High-Producer","Low-Producer","Non-Producer","Total Biomass","\beta-Carotene"],'Location','best');
    xlabel("Time (hr)"); ylabel("\beta-Carotene Concentration (mg/L)")
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

    prod_K_O2 = 0.05;
    
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
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(pi./24)*(cos(t./24.*pi));
    dydt(6) = p(21).*(y(5)./(prod_K_O2+y(5))).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        dV./V.*y(6);
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

    prod_K_O2 = 0.05;
    
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
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5)));
    % dydt(4) = 0;
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(cos(freq.*(2./pi).*t));
    dydt(6) = p(21).*(y(5)./(prod_K_O2+y(5))).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*(y(5)./(prod_K_O2+y(5))).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)));
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

    prod_K_O2 = 0.05;
    
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
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(cos(freq.*(2./pi).*t));
    dydt(6) = p(21).*(y(5)./(prod_K_O2+y(5))).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*(y(5)./(prod_K_O2+y(5))).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        dV./V.*y(6);
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

    prod_K_O2 = 0.05;

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
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(cos(freq.*(2./pi).*t));
    dydt(6) = p(21).*(y(5)./(prod_K_O2+y(5))).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*(y(5)./(prod_K_O2+y(5))).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        dV./V.*y(6);
    dydt(7) = dV.*y(6);

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
    
    prod_K_O2 = 0.05;
    
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
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(28).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) + dV./V.*C_G_in - dV./V.*y(4);
    % dydt(4) = 0;
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(29).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*(cos(freq.*(2./pi).*t)-y_shift);
    % dydt(5) = 0.007./4.*(pi./24)*(cos(t./24.*pi));
    dydt(6) = p(21).*(y(5)./(prod_K_O2+y(5))).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*(y(5)./(prod_K_O2+y(5))).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + ...
        p(30).*(y(5)./(prod_K_O2+y(5))).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - dV./V.*y(6);
    dydt(7) = dV.*y(6);
    dydt(8) = p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - p(26).*y(8) - ...
        p(27).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - dV./V.*y(8);
    % dydt(8) = 0;

    dydt = dydt';
end