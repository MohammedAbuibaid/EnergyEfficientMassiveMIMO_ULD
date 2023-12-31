%% Secnario A Results
clear all; close all; clc;
load('Res_Scenario_B.mat');

% Average EE Comparison 
figure
hold on; box on;
plot(1:24,AvgEEperBS/1e6,'.--g')
plot(1:24,AvgEEperBS_Fix/1e6,'.--r')
grid on
xlim([1 24])
ylim([0 18])
legend('ULD Adap. Sys', 'Fixed Anten. Sys','location','best')
ylabel('Energy Efficiency [Mbit/Joule]')
xlabel('Time: Hours');

% Average User Rate Comparison 
figure
hold on; box on; grid on;
set(gca,'fontsize',10)
plot(1:24,AvgURperBS/1e6,'.--g')
plot(1:24,AvgURperBS_Fix/1e6,'.--r')
xlim([1 24])
ylim([0 250])
legend('ULD Adap. Sys', 'Fixed Anten. Sys','location','best')
ylabel('User Rate [Mbps]')
xlabel('Time: Hours');

% Average Total Power Comparison 
figure
hold on; box on; grid on;
set(gca,'fontsize',10)
plot(1:24,AvgPtotperBS/1e3,'.--g')
plot(1:24,AvgPtotperBS_Fix/1e3,'.--r')
xlim([1 24])
ylim([0 2])
legend('ULD Adap. Sys', 'Fixed Anten. Sys','location','best')
ylabel('Average Total power [kW]')
xlabel('Time: Hours');

%%

% Hourly Num of Antennas
figure
hold on; box on; grid on;
set(gca,'fontsize',10)
plot(1:24,Mmax_h,'.--g')
plot(1:24,MgOpt*ones(1,24),'.--r')
xlim([1 24])
ylim([0 250])
legend('ULD Adap. Sys', 'Fixed Anten. Sys','location','best')
ylabel('Num. of Anten. [M]')
xlabel('Time: Hours');


%% Weighted Average Hourly Num of Antennas

figure
hold on; box on; grid on;
set(gca,'fontsize',10)
plot(1:24,WA_Mopt_h,'.--g')
plot(1:24,MgOpt*ones(1,24),'.--r')
xlim([1 24])
ylim([0 250])
legend('ULD Adap. Sys', 'Fixed Anten. Sys','location','best')
ylabel('WA Num. of Anten. [M]')
xlabel('Time: Hours');



%%





%% % Paper Figures 

load('Res_Scenario_A.mat');
Rmax = 0.5; Rmin = 0.035;
Cell_Area = 3*sqrt(3)*(Rmax)^2/2 - pi*(Rmin)^2;
Energy_per_Area_Fix = (mean(AvgPtotperBS_Fix)*60*60*(1/1e6))/Cell_Area;
Energy_per_Area_Ada = (mean(AvgPtotperBS)*60*60*(1/1e6))/Cell_Area;
disp(['Energy_per_Area_Fix  = ' num2str(Energy_per_Area_Fix) ' MWh/km2']);
disp(['Energy_per_Area_Ada  = ' num2str(Energy_per_Area_Ada) ' MWh/km2']);


%%
% Hourly Num of Antennas
load('Res_Scenario_A.mat');
Mavg_DH = MgOpt-mean(Mmax_h(8:21));
Mavg_NH = MgOpt - mean(Mmax_h([1:7 22:24]));
disp ('--------------------------------------------')
disp(['MgOpt = ' num2str(MgOpt)]);
disp ('For Scenario (A), Mavg can be turned OFF at:')
disp (['Daylight hours = ' num2str(Mavg_DH) ', Night hours = '  num2str(Mavg_NH)])

load('Res_Scenario_B.mat');
disp ('--------------------------------------------')
disp(['MgOpt = ' num2str(MgOpt)]);
disp ('For Scenario (B), Mavg can be turned OFF at:')
Mavg_DH = MgOpt-mean(Mmax_h(8:21));
Mavg_NH = MgOpt - mean(Mmax_h([1:7 22:24]));
disp (['Daylight hours = ' num2str(Mavg_DH) ', Night hours = '  num2str(Mavg_NH)])
disp ('--------------------------------------------')

%% 
figure
hold on; box on; 
set(gca,'fontsize',12)

load('Res_Scenario_B.mat');
plot(1:24,MgOpt*ones(1,24),'-r','linewidth',1.2)

load('Res_Scenario_A.mat');
plot(1:24,MgOpt*ones(1,24),'-b','linewidth',1.2)


load('Res_Scenario_B.mat');
plot(1:24,Mmax_h,'--or','linewidth',1.2 )

load('Res_Scenario_A.mat');
plot(1:24,Mmax_h,'--xb','linewidth',1.2 )

xlim([1 24])
set(gca, 'XTick', 2:2:24)
ylim([50 250])

legend('Fixed Sys. Scenario B','Fixed Sys. Scenario A','Adap. Sys. Scenario B','Adap. Sys. Scenario A' ,'location', 'best')

ylabel('Max Number of Antennas')
xlabel('Time: Hours');


%% AVG M


figure
hold on; box on; 
set(gca,'fontsize',12)

load('Res_Scenario_B.mat');
plot(1:24,MgOpt*ones(1,24),'-r','linewidth',1.2)

load('Res_Scenario_A.mat');
plot(1:24,MgOpt*ones(1,24),'-b','linewidth',1.2)


load('Res_Scenario_B.mat');
plot(1:24,WA_Mopt_h,'--or','linewidth',1.2 )

load('Res_Scenario_A.mat');
plot(1:24,WA_Mopt_h,'--xb','linewidth',1.2 )

xlim([1 24])
set(gca, 'XTick', 2:2:24)
ylim([0 250])

legend('Fixed Sys. Scenario B','Fixed Sys. Scenario A','Adap. Sys. Scenario B','Adap. Sys. Scenario A' ,'location', 'best')

ylabel('Average Number of Antennas')
xlabel('Time: Hours');

set(gcf, 'Position', [500 500 900 500])




%% Hourly Num of Antennas
load('Res_Scenario_A.mat');
Mavg_DH = MgOpt - mean(WA_Mopt_h(8:21));
Mavg_NH = MgOpt - mean(WA_Mopt_h([1:7 22:24]));
disp ('--------------------------------------------')
disp(['MgOpt = ' num2str(MgOpt)]);
disp ('For Scenario (A), Mavg can be turned OFF at:')
disp (['Daylight hours = ' num2str(Mavg_DH) ', Night hours = '  num2str(Mavg_NH)])

load('Res_Scenario_B.mat');
disp ('--------------------------------------------')
disp(['MgOpt = ' num2str(MgOpt)]);
disp ('For Scenario (B), Mavg can be turned OFF at:')
Mavg_DH = MgOpt - mean(WA_Mopt_h(8:21));
Mavg_NH = MgOpt - mean(WA_Mopt_h([1:7 22:24]));
disp (['Daylight hours = ' num2str(Mavg_DH) ', Night hours = '  num2str(Mavg_NH)])
disp ('--------------------------------------------')



%% Hourly Num of Antennas
load('Res_Scenario_A.mat');
Mavg_24hr = MgOpt - mean(WA_Mopt_h);
disp ('--------------------------------------------')
disp(['MgOpt = ' num2str(MgOpt)]);
disp ('For Scenario (A), Mavg can be turned OFF at:')
disp (['24 hours = ' num2str(Mavg_24hr)])
Mavg_24hr/MgOpt


load('Res_Scenario_B.mat');
disp ('--------------------------------------------')
disp(['MgOpt = ' num2str(MgOpt)]);
disp ('For Scenario (B), Mavg can be turned OFF at:')
Mavg_24hr = MgOpt - mean(WA_Mopt_h);
disp (['24 hours = ' num2str(Mavg_24hr)])
disp ('--------------------------------------------')
Mavg_24hr/MgOpt


