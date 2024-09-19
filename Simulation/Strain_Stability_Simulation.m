clear;
clc;
close all;

%% Setup
mu_max_1 = 0.3466.*0.6;
K_HI_G = 1;
K_HI_O2 = 0.04;
k_d_HI = 0.001;
MR_HI = 0.1;

mu_max_2 = 0.3466.*0.8;
K_LO_G = 1;
K_LO_O2 = 0.04;
k_d_LO = 0.001;
MR_LO = 0.1;

mu_max_3 = 0.3466;
K_NO_G = 1;
K_NO_O2 = 0.04;
k_d_NO = 0.001;

Y_HI_G = 2;
Y_LO_G = 2;
Y_NO_G = 2;

Y_HI_O2 = 0.05;
Y_LO_O2 = 0.05;
Y_NO_O2 = 0.05;

Y_HI_prod = 0.1;
Y_LO_prod = 0.05;

mu_max_4 = 0.3466;
K_NOF_G = 1;
K_NOF_O2 = 0.04;
k_d_NOF = 0.001;

p = [mu_max_1,K_HI_G,K_HI_O2,k_d_HI,MR_HI,mu_max_2,K_LO_G,K_LO_O2,k_d_LO,MR_LO, ...
        mu_max_3,K_NO_G,K_NO_O2,k_d_NO,Y_HI_G,Y_LO_G,Y_NO_G,Y_HI_O2,Y_LO_O2,Y_NO_O2, ...
        Y_HI_prod,Y_LO_prod,mu_max_4,K_NOF_G,K_NOF_O2,k_d_NOF];

%% Run Model
batch_ct = 1;
t_res = [];
y_res = cell(6,1);
tot_prod_mass = 0;
V = 100; % L
batch_time_length = 50; % hr
t_span_1 = [0,batch_time_length];
ODE_1 = ode;
ODE_1.ODEFcn = @(t,y,p) ODESys_batch(t,y,p);
ODE_1.InitialTime = 0;
ODE_1.InitialValue = [0.1,0.025,0,10,0.0416,0];
ODE_1.Parameters = p;
ODE_1.Solver = 'ode15s';

sol = solve(ODE_1,t_span_1(1),t_span_1(2));
t_res = sol.Time;
y_res = sol.Solution;

fig1 = figure;

subplot(2,1,1);
hold on;
plot(t_res,y_res(1,:),'b-','LineWidth',2);
plot(t_res,y_res(2,:),'k-','LineWidth',2);
plot(t_res,y_res(3,:),'r-','LineWidth',2);
legend(["High-Producer","Low-Producer","Non-Producer"]);
xlabel("Time (hr)"); ylabel("Biomass Concentration (g/L)");
hold off;

subplot(2,1,2);
hold on;
yyaxis left;
plot(t_res,y_res(4,:),'b-','LineWidth',2);
yyaxis right;
plot(t_res,y_res(6,:),'k-','LineWidth',2);
legend("Product");
xlabel("Time (hr)"); ylabel("Concentration (g/L)");
hold off;

axes = fig1.Children;
for k=1:1:length(axes)
    axes(k).FontName = 'Arial';
    axes(k).FontSize= 12;
    axes(k).FontWeight = 'bold';
end

%% Functions
function dydt = ODESys_batch(t,y,p)
    % 1. HI
    % 2. LO
    % 3. NO
    % 4. Glucose
    % 5. O2
    % 6. prod
    % 7. prod mass

    S_O2 = 0.0416;
    k_L_a = 100;
    amp = 0.03;
    freq = 0.1;

    dydt(1) = p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(4).*y(1) - p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)));
    dydt(2) = p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - p(9).*y(2) - ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)));
    dydt(3) = p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) - p(14).*y(3) + ...
        p(10).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) + p(5).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5)));
    dydt(4) = -p(15).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(16).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(17).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5)));
    dydt(5) = -p(18).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) - p(19).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5))) - ...
        p(20).*p(11).*y(3).*(y(4)./(p(12)+y(4))).*(y(5)./(p(13)+y(5))) + k_L_a.*(S_O2-y(5)) + amp.*cos(freq.*(2./pi).*t);
    dydt(6) = p(21).*p(1).*y(1).*(y(4)./(p(2)+y(4))).*(y(5)./(p(3)+y(5))) + p(22).*p(6).*y(2).*(y(4)./(p(7)+y(4))).*(y(5)./(p(8)+y(5)));

    dydt = dydt';
end