% Bladder growth model framework
% --------------- WORKSPACE CLEANUP -------------
clear
clc
close all
global k_collagen k_elastin radius_tzero thickness_tzero ratio deltaT% Global variable  k_collagen(collagen stiffness)

% -------------------------------------------------------------------------
%% Bladder dimensions.
%load the recruitment stretch: min mod max
load recruitment_stretch.mat
R_LP{1} = LR{4};
R_DSM{1} = DR{4};
%load the wall thickness t and radius r.
r = [6.48 6.87 6.82 7.03 6.27 7.43 7.76 6.60 6.04];
t = [0.66 0.52 0.61 0.46 0.75 0.85 0.92 1.12 0.66];
% Unloaded radius at t = 0 (m)
radius_tzero    = r(4)*1e-3;
% Unloaded thickness at t=0 (m).
thickness_tzero = t(4)*1e-3;
%% Define material constants
% Define elastin stiffness
k_elastin = 63.3055;%375.9471;
% Define collagen stiffness
k_collagen = 4.3755e5;%858780;
% Collagen ratio between lamina priopria and detrusor
ratio = 4;
% Threshold pressure (Pa) (=15cmH2O)
pressure = 2500;
stress = 100e3;
% Define the thick_factor
thick_factor(1) = 1;
Mass(1) =1;
% Remodeling loop
n = 100;
for i = 1:n
% if n<21
%     pressure = 2500;
% else
%     pressure = 5000;
% end

%% Calculate threshold stretch
lambda_th(i) = solv(pressure,R_LP{i},R_DSM{i});
%% Define collagen and smooth muscle cells attachment distribution
% Lamina propria (LP) collagen attachment stretch distribution
[AT_LP{1},AT_LP_W(1),AT_LP_S(1)] = collagen_stretch(R_LP{1},lambda_th(1));
% Destrusor (DSM) collagen attachment stretch
[AT_DSM{1},AT_DSM_W(1),AT_DSM_S(1)] = collagen_stretch(R_DSM{1},lambda_th(1));
% Evolution of collagen attachment distribution

if 1<i<=10
AT_DSM{i+1}(1) = AT_DSM{i}(1);
AT_DSM{i+1}(3) = AT_DSM{i}(3);
AT_DSM{i+1}(2) = AT_DSM{i}(2);
else
AT_DSM{i+1}(1) = AT_DSM{i}(1);
AT_DSM{i+1}(3) = AT_DSM{i}(3);
AT_DSM{i+1}(2) = AT_DSM{i}(2);    
end

%% ---------------- Calculate stress & pressure at very time point ---------------#
% Set up the range of stretch
lam{i} = 1:0.001:lambda_th(i);
% Calculate stress
[S_C_LP{i}, S_C_DSM{i}, S_E{i}, S{i}] = stress_cal(R_DSM{i},R_LP{i},k_collagen,k_elastin,lam{i});
% Define dimension term
D =  2.*thickness_tzero ./ (radius_tzero .* lam{i}.^3);
% Pressure on LP collagen
P_C_LP{i} = S_C_LP{i}.*D;
% Pressure on DSM collagen
P_C_DSM{i} = S_C_DSM{i}.*D;
% Pressue on elastin;
P_E{i} = S_E{i}.*D;
% Overall pressure
P{i} = S{i}.*D;

%% ---------------- Collagen remodeling ---------------#
% Initial mass 
M_LP(1) = 1;
M_DSM(1) = 1;
% Collagen remodeling & growth rate
Alpha_LP = 0.2;%*(max(P_C_LP{i})-max(P_C_LP{1}))/(4000); % Remodeling constant
Beta_LP = 0.2; % Growth constant
Alpha_DSM = 0.2;%*(max(P_C_DSM{i})-max(P_C_DSM{1}))/(4000); % Remodeling constant
Beta_DSM = 0.6; % Growth constant
% Collagen remodeling time step
deltaT = 1;

%% SMC death change collagen distribution 

% if i == 2
% R_DSM{i}(1) = R_DSM{i}(1)-0.2;%-0.01
% R_DSM{i}(3) = R_DSM{i}(3);%-0.01
% R_DSM{i}(2) = R_DSM{i}(2);
% end


% pressure = pressure+20;



% Lamina propria collagen remodeling
[R_LP{i+1},M_LP(i+1)] = collagen_remodeling(M_LP(i),R_LP{i},AT_LP{1},Alpha_LP,Beta_LP,lambda_th(i));
% Detrusor collagen remodeling
[R_DSM{i+1},M_DSM(i+1)] = collagen_remodeling(M_DSM(i),R_DSM{i},AT_DSM{1},Alpha_DSM,Beta_DSM,lambda_th(i));


% %% Growth
% thick_factor(i+1) = 0.5*M_DSM(i+1)+0.5*M_LP(i+1);
% 
% thick_factor(i+1) = thick_factor(i);
% 
% Mass(i+1) = 0.5*(M_DSM(i+1))+0.5*M_LP(i+1);
% 
% k_collagen = 4.3755e5*M_DSM(i);
%  
% ratio = 4*M_LP(i+1)/M_DSM(i+1);
% thickness_tzero = thickness_tzero*1.03; 

end

%% ---------------- Plot pdf/cdf distributon ---------------#
collagen_disp(R_LP{1},R_DSM{1},1,2);
collagen_disp(R_LP{n},R_DSM{n},3,4);

%% ---------------- Plot the loading curve (Mechanical properties) ---------------#
% Plot loading curve of intial and end point
figure(5)
hold on
plot(lam{1},S{1},'LineWidth',3)
plot(lam{n},S{n},'LineWidth',3)
legend('t=1','t=n')
xlabel('Stretch')
ylabel('Stress (KPa)')
set(gca,'fontsize',15)
grid
%% ---------------- Plot the Pressure inflation curve ---------------#
% Plot the initial time point
figure(6)
hold on
plot(lam{1},P{1},'LineWidth',3)
plot(lam{1},P_C_LP{1},'LineWidth',3)
plot(lam{1},P_C_DSM{1},'LineWidth',3)
plot(lam{1},P_E{1},'LineWidth',3)
legend('Overall pressure', 'Lamina propria pressure', 'Detrusor pressure', 'Elastin pressure')
xlabel('Stretch')
ylabel('Pressure (KPa)')
set(gca,'fontsize',15)
grid
% Plot the end time point
figure(7)
hold on
plot(lam{n},P{n},'LineWidth',3)
plot(lam{n},P_C_LP{n},'LineWidth',3)
plot(lam{n},P_C_DSM{n},'LineWidth',3)
plot(lam{n},P_E{n},'LineWidth',3)
legend('Overall pressure', 'Lamina propria pressure', 'Detrusor pressure', 'Elastin pressure')
xlabel('Stretch')
ylabel('Pressure (KPa)')
set(gca,'fontsize',15)
grid
%  %% ---------------- Plot the mass change ---------------#
% figure(8)
% hold on
% plot(1:n+1,M_LP,'LineWidth',3)
% plot(1:n+1,M_DSM,'LineWidth',3)
% plot(1:n+1,Mass,'LineWidth',3)
% legend('LP collagen mass', 'DS collagen mass','Total mass')
% xlabel('Time')
% ylabel('Mass density')
% set(gca,'fontsize',15)
% grid
%  %% ---------------- Plot the thickness ---------------#
% figure(9)
% hold on
% plot(1:n+1,thickness_tzero*thick_factor,'LineWidth',3)
% xlabel('Time')
% ylabel('thickness')
% set(gca,'fontsize',15)
% grid
