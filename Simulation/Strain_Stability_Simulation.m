clear;
clc;
close all;

%% Tasks:
% 1. Higher dilution rate => faster loss of \beta-Caroteneion
% 2. Increased heterogeneity => slower loss of \beta-Caroteneion
% 3. Oxygen concentration impact on growth rate?

%% Setup
mu_max_1 = 0.3466.*0.6;
K_HI_G = 1;
K_HI_O2 = 0.04;
k_d_HI = 0.001;
MR_HI = 0.001;

mu_max_2 = 0.3466.*0.95;
K_LO_G = 1;
K_LO_O2 = 0.04;
k_d_LO = 0.001;
MR_LO = 0.001;

mu_max_3 = 0.3466;
K_NO_G = 1;
K_NO_O2 = 0.04;
k_d_NO = 0.001;

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

mu_max_4 = 0.3466.*0.9;
K_MED_G = 1;
K_MED_O2 = 0.04;
k_d_MED = 0.001;
MR_MED = 0.001;

p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
        mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
        Y_HI_prod,Y_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,Y_MED_prod];

param_names = "mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO,mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2,Y_HI_prod,Y_LO_prod,mu_max_4,K_MED_G,K_MED_O2,k_d_MED,MR_MED,Y_MED_G,Y_MED_O2,Y_MED_prod";
param_names = split(param_names,',');

disp(table((1:1:length(p))',p',param_names,'VariableNames',{'Number','Value','Name'}))

% Environmental Conditions
dV = 5; % L/s
V = 100;% L
C_G_in = 10; % g/L
S_O2 = 0.0416;
k_L_a = 1;
amp = 0;
freq = 0.1;
env_cond = [dV,V,C_G_in,S_O2,k_L_a,amp,freq];

%% Run Model
% Task 1
fig1 = figure;
fig1.Color = [1,1,1];

cultivation_time = 1000; % hr
t_span = [0,cultivation_time];

% base dilution rate (0.05 1/hr)
ODE_1 = ode;
ODE_1.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
ODE_1.InitialTime = 0;
ODE_1.InitialValue = [1,0,0,10,0.0416,0,0];
ODE_1.Parameters = p;
ODE_1.Solver = 'ode15s';

sol = solve(ODE_1,t_span(1),t_span(2));
t_res = sol.Time;
y_res = sol.Solution;

subplot(3,2,1);
hold on;
plot(t_res,y_res(1,:),'b-','LineWidth',2);
plot(t_res,y_res(2,:),'k-','LineWidth',2);
plot(t_res,y_res(3,:),'r-','LineWidth',2);
legend(["High-Producer","Low-Producer","Non-Producer"],'Location','east');
xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
title("Dilution Rate: 0.05 1/{hr}");
hold off;

subplot(3,2,3);
hold on;
yyaxis left;
plot(t_res,y_res(4,:),'b-','LineWidth',2);
ylabel("Glucose");
axes = fig1.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end
0;
plot(t_res,y_res(6,:),'k-','LineWidth',2);
% plot(t_res,y_res(7,:),'r-','LineWidth',2);
ylim([0,0.5]);
legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (g/L)"],'Location','east');
xlabel("Time (hr)"); ylabel("\beta-Carotene");
hold off;

subplot(3,2,5);
hold on;
plot(t_res,y_res(5,:),'b-','LineWidth',2);
xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
ylim([0,Inf]);
hold off;

% higher dilution rate (0.1 1/hr)
env_cond(1) = 10; % 1/hr
ODE_2 = ode;
ODE_2.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
ODE_2.InitialTime = 0;
ODE_2.InitialValue = [1,0.25,0,10,0.0416,0,0];
ODE_2.Parameters = p;
ODE_2.Solver = 'ode15s';
env_cond(1) = 5; % 1/hr

sol = solve(ODE_2,t_span(1),t_span(2));
t_res_2 = sol.Time;
y_res_2 = sol.Solution;

subplot(3,2,2);
hold on;
plot(t_res_2,y_res_2(1,:),'b-','LineWidth',2);
plot(t_res_2,y_res_2(2,:),'k-','LineWidth',2);
plot(t_res_2,y_res_2(3,:),'r-','LineWidth',2);
legend(["High-Producer","Low-Producer","Non-Producer"],'Location','east');
xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
title("Dilution Rate: 0.1 1/{hr}");
ylim([0,5]);
hold off;

subplot(3,2,4);
hold on;
yyaxis left;
plot(t_res_2,y_res_2(4,:),'b-','LineWidth',2);
ylabel("Glucose"); ylim([0,Inf]);
axes = fig1.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end
yyaxis right;
plot(t_res_2,y_res_2(6,:),'k-','LineWidth',2);
% plot(t_res_2,y_res_2(7,:),'r-','LineWidth',2);
ylim([0,0.5]);
legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (g/L)"],'Location','east');
xlabel("Time (hr)"); ylabel("\beta-Carotene");
hold off;

subplot(3,2,6);
hold on;
plot(t_res_2,y_res_2(5,:),'b-','LineWidth',2);
xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
ylim([0,Inf]);
hold off;

axes = fig1.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end

% Task 2
fig2 = figure;
fig2.Color = [1,1,1];

cultivation_time = 1000; % hr
t_span = [0,cultivation_time];

% 2 different strains
ODE_3 = ode;
ODE_3.ODEFcn = @(t,y,p) ODESys_cont_2(t,y,p,env_cond);
ODE_3.InitialTime = 0;
ODE_3.InitialValue = [1.5,0,0,10,0.0416,0,0];
ODE_3.Parameters = p;
ODE_3.Solver = 'ode15s';

sol = solve(ODE_3,t_span(1),t_span(2));
t_res = sol.Time;
y_res = sol.Solution;

subplot(3,3,1);
hold on;
plot(t_res,y_res(1,:),'b-','LineWidth',2);
plot(t_res,y_res(3,:),'r-','LineWidth',2);
legend(["High-Producer","Non-Producer"],'Location','east');
xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
title("2 Distinct Strains");
ylim([0,5]);
hold off;

subplot(3,3,4);
hold on;
yyaxis left;
plot(t_res,y_res(4,:),'b-','LineWidth',2);
ylabel("Glucose");
axes = fig2.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end
yyaxis right;
plot(t_res,y_res(6,:),'k-','LineWidth',2);
% plot(t_res,y_res(7,:),'r-','LineWidth',2);
ylim([0,0.5]);
legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (g/L)"],'Location','east');
xlabel("Time (hr)"); ylabel("\beta-Carotene");
hold off;

subplot(3,3,7);
hold on;
plot(t_res,y_res(5,:),'b-','LineWidth',2);
xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
ylim([0,Inf]);
hold off;

% 3 different strains
ODE_4 = ode;
ODE_4.ODEFcn = @(t,y,p) ODESys_cont_3(t,y,p,env_cond);
ODE_4.InitialTime = 0;
ODE_4.InitialValue = [1.25,0.25,0,10,0.0416,0,0];
ODE_4.Parameters = p;
ODE_4.Solver = 'ode15s';

sol = solve(ODE_4,t_span(1),t_span(2));
t_res_2 = sol.Time;
y_res_2 = sol.Solution;

subplot(3,3,2);
hold on;
plot(t_res_2,y_res_2(1,:),'b-','LineWidth',2);
plot(t_res_2,y_res_2(2,:),'k-','LineWidth',2);
plot(t_res_2,y_res_2(3,:),'r-','LineWidth',2);
legend(["High-Producer","Low-Producer","Non-Producer"],'Location','east');
xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
title("3 Distinct Strains");
ylim([0,5]);
hold off;

subplot(3,3,5);
hold on;
yyaxis left;
plot(t_res_2,y_res_2(4,:),'b-','LineWidth',2);
ylabel("Glucose"); ylim([0,Inf]);
axes = fig2.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end
yyaxis right;
plot(t_res_2,y_res_2(6,:),'k-','LineWidth',2);
% plot(t_res_2,y_res_2(7,:),'r-','LineWidth',2);
ylim([0,0.5]);
legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (g/L)"],'Location','east');
xlabel("Time (hr)"); ylabel("\beta-Carotene");
hold off;

subplot(3,3,8);
hold on;
plot(t_res_2,y_res_2(5,:),'b-','LineWidth',2);
xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
ylim([0,Inf]);
hold off;

% 4 different strains
ODE_5 = ode;
ODE_5.ODEFcn = @(t,y,p) ODESys_cont_4(t,y,p,env_cond);
ODE_5.InitialTime = 0;
ODE_5.InitialValue = [1,0.25,0.25,10,0.0416,0,0,0];
ODE_5.Parameters = p;
ODE_5.Solver = 'ode15s';

sol = solve(ODE_5,t_span(1),t_span(2));
t_res_3 = sol.Time;
y_res_3 = sol.Solution;

subplot(3,3,3);
hold on;
plot(t_res_3,y_res_3(1,:),'b-','LineWidth',2);
plot(t_res_3,y_res_3(8,:),'g-','LineWidth',2);
plot(t_res_3,y_res_3(2,:),'k-','LineWidth',2);
plot(t_res_3,y_res_3(3,:),'r-','LineWidth',2);
legend(["High-Producer","Medium-Producer","Low-Producer","Non-Producer"],'Location','east');
xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
title("4 Distinct Strains");
ylim([0,5]);
hold off;

subplot(3,3,6);
hold on;
yyaxis left;
plot(t_res_3,y_res_3(4,:),'b-','LineWidth',2);
ylabel("Glucose"); ylim([0,Inf]);
axes = fig2.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end
yyaxis right;
plot(t_res_3,y_res_3(6,:),'k-','LineWidth',2);
% plot(t_res_3,y_res_3(7,:),'r-','LineWidth',2);
ylim([0,0.5]);
legend(["Glucose Concentration (g/L)","\beta-Carotene Concentration (g/L)"],'Location','east');
xlabel("Time (hr)"); ylabel("\beta-Carotene");
hold off;

subplot(3,3,9);
hold on;
plot(t_res_3,y_res_3(5,:),'b-','LineWidth',2);
xlabel("Time (hr)"); ylabel("DO Concentration (g/L)");
ylim([0,Inf]);
hold off;

axes = fig2.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end

disp(table(y_res(7,end),y_res_2(7,end),y_res_3(7,end),'VariableNames',{'2','3','4'}))

% Task 3
fig3 = figure;
fig3.Color = [1,1,1];
ODE_6 = ode;
ODE_6.ODEFcn = @(t,y,p) simple_mu_exper(t,y,p);
ODE_6.InitialTime = 0;
ODE_6.InitialValue = [1,0.01,10];
params = [[1,2].*0.08,[1,1,0.5,0.5]];
ODE_6.Parameters = params;
ODE_6.Solver = 'ode15s';

t_span = [0,500];
sol = solve(ODE_6,t_span(1),t_span(2));
t_res = sol.Time;
y_res = sol.Solution;

subplot(1,3,1);
hold on;
plot(t_res,ones(size(t_res)).*params(1),'b-','LineWidth',2);
plot(t_res,ones(size(t_res)).*params(2),'k-','LineWidth',2);
xlabel("Time (hr)"); ylabel("Specific Maximum Growth Rate (1/hr)");
legend(["ES","WT"],'Location','northeast');
ylim([0,0.3]);
hold off;

subplot(1,3,2);
hold on;
plot(t_res,params(1).*y_res(1,:).*(y_res(3,:)./(p(3)+y_res(3,:))),'b-','LineWidth',2);
plot(t_res,params(2).*y_res(2,:).*(y_res(3,:)./(p(4)+y_res(3,:))),'k-','LineWidth',2);
xlabel("Time (hr)"); ylabel("Growth Rate (g/(L*hr))");
legend(["ES","WT"],'Location','east');
hold off;

subplot(1,3,3);
hold on;
plot(t_res,y_res(1,:),'b-','LineWidth',2);
plot(t_res,y_res(2,:),'k-','LineWidth',2);
xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
legend(["ES","WT"],'Location','east');
hold off;

axes = fig3.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end

% closer mu
fig4 = figure;
fig4.Color = [1,1,1];
ODE_6 = ode;
ODE_6.ODEFcn = @(t,y,p) simple_mu_exper(t,y,p);
ODE_6.InitialTime = 0;
ODE_6.InitialValue = [1,0.01,10];
params = [[1.5,2].*0.08,[1,1,0.5,0.5]];
ODE_6.Parameters = params;
ODE_6.Solver = 'ode15s';

sol = solve(ODE_6,t_span(1),t_span(2));
t_res = sol.Time;
y_res = sol.Solution;

subplot(1,3,1);
hold on;
plot(t_res,ones(size(t_res)).*params(1),'b-','LineWidth',2);
plot(t_res,ones(size(t_res)).*params(2),'k-','LineWidth',2);
xlabel("Time (hr)"); ylabel("Specific Maximum Growth Rate (1/hr)");
legend(["ES","WT"],'Location','northeast');
ylim([0,0.3]);
hold off;

subplot(1,3,2);
hold on;
plot(t_res,params(1).*y_res(1,:).*(y_res(3,:)./(p(3)+y_res(3,:))),'b-','LineWidth',2);
plot(t_res,params(2).*y_res(2,:).*(y_res(3,:)./(p(4)+y_res(3,:))),'k-','LineWidth',2);
xlabel("Time (hr)"); ylabel("Growth Rate (g/(L*hr))");
legend(["ES","WT"],'Location','east');
ylim([0,3.5]);
hold off;

subplot(1,3,3);
hold on;
plot(t_res,y_res(1,:),'b-','LineWidth',2);
plot(t_res,y_res(2,:),'k-','LineWidth',2);
xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
legend(["ES","WT"],'Location','east');
hold off;

axes = fig4.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 8;
    axes(k).FontWeight = 'bold';
    axes(k).Box = true;
    try
        axes(k).YColor = [0,0,0];
    catch err

    end
end

%% Functions
function dydt = simple_mu_exper(t,y,p)
    dV = 5;
    V = 100;
    C_in = 10;
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
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq] = deal(env_cond_cell{:});

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
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*cos(freq.*(2./pi).*t);
    % dydt(5) = 0.0416./4.*(pi./24)*(cos(t./24.*pi));
    dydt(6) = p(21).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
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

    env_cond_cell = num2cell(env_cond);
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq] = deal(env_cond_cell{:});

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
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*cos(freq.*(2./pi).*t);
    % dydt(5) = 0.0416./4.*(pi./24)*(cos(t./24.*pi));
    dydt(6) = p(21).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
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

    env_cond_cell = num2cell(env_cond);
    [dV,V,C_G_in,S_O2,k_L_a,amp,freq] = deal(env_cond_cell{:});

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
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(29).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*cos(freq.*(2./pi).*t);
    % dydt(5) = 0.0416./4.*(pi./24)*(cos(t./24.*pi));
    dydt(6) = p(21).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + ...
        p(30).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - dV./V.*y(6);
    dydt(7) = dV.*y(6);
    dydt(8) = p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - p(26).*y(8) - ...
        p(27).*p(23).*y(8).*(y(4)./(p(24)+y(4))).*(y(5)./(p(25)+y(5))) - dV./V.*y(8);
    % dydt(8) = 0;

    dydt = dydt';
end