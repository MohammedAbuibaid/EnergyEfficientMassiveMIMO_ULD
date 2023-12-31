clc; clear all; close all; 
% Simulation Enviroment
Rmax = 500; %Cell radius (distance to a vertex of hexagonal cell)
Rmin = 35; %Users inside this circle will not considered in simulations
TestPoints = 10000;
ISD = Rmax*sqrt(3);
% % Coordinates of all BSs in the Hexagonal Network
u = [0 1 0 -1 -1  0  1  2 2 1 0 -1 -2 -2 -2 -1  0  1  2]; % 30-Degree axis
v = [0 0 1  1  0 -1 -1 -1 0 1 2  2  2  1  0 -1 -2 -2 -2]; % Vertical axis
% Setting the BSs in their locations as seen from the origin.
BSLocations = sqrt(3).*(ISD/2+1i*Rmax/2).*u + (0+1i*ISD).*v;
Pc = 20;
Mmax = 250;Kmax=125;
PmaxPA = 10^0.6; % 6dB
p_max = PmaxPA/10^0.8;
Mmin = ceil(Pc/p_max);
GoS = 0.02;
loading = [10 20 30 40 50 60 70 80 90 100];




%% Performance of Adaptive System in BF
load('Res_BF.mat');
Ri = [Rmin 200 400 Rmax];
ULD_BF = [0.10 0.20 0.70]';
UELocations = UE_insertion_MonteCarlo_HexCell(TestPoints,ULD_BF,Ri,false);
[PLO_Network,PLI_Network] = Wrap_Around_PLO_PLI(BSLocations,UELocations,1,Rmax,false);

[BF_EE_Fix,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,MgOpt*ones(1,KgOpt),Pc,PmaxPA,false);
[BF_EE_Lo10,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,BF_Mopt_Lo(1,:),Pc,PmaxPA,false);
[BF_EE_Lo50,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,BF_Mopt_Lo(5,:),Pc,PmaxPA,false);
[BF_EE_Lo100,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,BF_Mopt_Lo(10,:),Pc,PmaxPA,false);

% figure % 1
% hold on; grid on; box on
% plot(BF_EE_Fix/1e6,'-.b')
% plot(BF_EE_Lo10/1e6,'.-g')
% plot(BF_EE_Lo50/1e6,'.-k')
% plot(BF_EE_Lo100/1e6,'.-r')
% xlim([0 110])
% title('Boundary Focused User Dist')
% legend('Fix MgOpt','Load 10%','Load 50%','Load 100%','location','best')
% ylabel('Energy Efficiency, [Mbit/Joule]');
% xlabel('Number of active users in the cell, K');



%% Performance of Adaptive System in Uniform ULD

load('Res_Uniform.mat');
Ri = [Rmin Rmax];
ULD_Uni = 1 ;
UELocations = UE_insertion_MonteCarlo_HexCell(TestPoints,ULD_Uni,Ri,false);
[PLO_Network,PLI_Network] = Wrap_Around_PLO_PLI(BSLocations,UELocations,1,Rmax,false);

[Uni_EE_Fix,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,MgOpt*ones(1,KgOpt),Pc,PmaxPA,false);
[Uni_EE_Lo10,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,Uni_Mopt_Lo(1,:),Pc,PmaxPA,false);
[Uni_EE_Lo50,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,Uni_Mopt_Lo(5,:),Pc,PmaxPA,false);
[Uni_EE_Lo100,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,Uni_Mopt_Lo(10,:),Pc,PmaxPA,false);

% figure %2
% hold on; grid on; box on
% plot(Uni_EE_Fix/1e6,'-.b')
% plot(Uni_EE_Lo10/1e6,'.-g')
% plot(Uni_EE_Lo50/1e6,'.-k')
% plot(Uni_EE_Lo100/1e6,'.-r')
% xlim([0 110])
% title('Uniform User Dist')
% legend('Fix MgOpt','Load 10%','Load 50%','Load 100%','location','best')
% ylabel('Energy Efficiency, [Mbit/Joule]');
% xlabel('Number of active users in the cell, K');



%% Performance of Adaptive System in CF ULD
load('Res_CF.mat');
Ri = [Rmin 200 400 Rmax];
ULD_CF = [0.80 0.10 0.10]';
UELocations = UE_insertion_MonteCarlo_HexCell(TestPoints,ULD_CF,Ri,false);
[PLO_Network,PLI_Network] = Wrap_Around_PLO_PLI(BSLocations,UELocations,1,Rmax,false);

[CF_EE_Fix,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,MgOpt*ones(1,KgOpt),Pc,PmaxPA,false);
[CF_EE_Lo10,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,CF_Mopt_Lo(1,:),Pc,PmaxPA,false);
[CF_EE_Lo50,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,CF_Mopt_Lo(5,:),Pc,PmaxPA,false);
[CF_EE_Lo100,~] = EE_R_for_Mopt(PLO_Network,PLI_Network,KgOpt,CF_Mopt_Lo(10,:),Pc,PmaxPA,false);

% figure %3
% hold on; grid on; box on
% plot(CF_EE_Fix/1e6,'-.b')
% plot(CF_EE_Lo10/1e6,'.-g')
% plot(CF_EE_Lo50/1e6,'.-k')
% plot(CF_EE_Lo100/1e6,'.-r')
% xlim([0 110])
% title('Center Focused User Dist')
% legend('Fix MgOpt','Load 10%','Load 50%','Load 100%','location','best')
% ylabel('Energy Efficiency, [Mbit/Joule]');
% xlabel('Number of active users in the cell, K');
 

%%

figure (4) 

hold on; grid on; box on
plot(BF_EE_Fix/1e6,'-.b')
plot(BF_EE_Lo10/1e6,'.-g')
plot(BF_EE_Lo50/1e6,'.-k')
plot(BF_EE_Lo100/1e6,'.-r')

plot(Uni_EE_Fix/1e6,'-.b')
plot(Uni_EE_Lo10/1e6,'.-g')
plot(Uni_EE_Lo50/1e6,'.-k')
plot(Uni_EE_Lo100/1e6,'.-r')

plot(CF_EE_Fix/1e6,'-.b')
plot(CF_EE_Lo10/1e6,'.-g')
plot(CF_EE_Lo50/1e6,'.-k')
plot(CF_EE_Lo100/1e6,'.-r')

xlim([0 110])
legend('Fix MgOpt','Load 10%','Load 50%','Load 100%','location','best')
ylabel('Energy Efficiency, [Mbit/Joule]');
xlabel('Number of active users in the cell, K');


%% BF Adaptive System: Mopt, EE and UR 

% figure %5
% hold on; grid on; box on
% plot(MgOpt*ones(1,KgOpt),'-.b')
% plot(BF_Mopt_Lo(1,:),'.-g')
% plot(BF_Mopt_Lo(5,:),'.-k')
% plot(BF_Mopt_Lo(10,:),'.-r')
% xlim([0 110])
% ylim([20 250])
% title('Boundary Focused User Dist')
% legend('Fix MgOpt','Load 10%','Load 50%','Load 100%','location','best')
% ylabel('Number of Active Antennas, M');
% xlabel('Number of active users in the cell, K');

figure%6
hold on; grid on; box on
plot(loading,BF_AvgEEperBS/1e6,'-xg')
plot(loading,BF_AvgEEperBS_Fix/1e6,'-xr')
ylim([0 20])
xlim([0 100])
title('Boundary Focused User Dist')
ylabel('Energy Efficiency [Mbits/Joule]');
xlabel('Network load (%)')
legend('Adaptive','Fix MgOpt','location','best')


figure %7
hold on; grid on; box on
plot(loading,BF_AvgURperBS/1e6,'-xg')
plot(loading,BF_AvgURperBS_Fix/1e6,'-xr')
ylim([0 250])
xlim([0 100])
title('Boundary Focused User Dist')
ylabel('Average rate per User [Mbps]');
xlabel('Network load (%)');
legend('Adaptive','Fix MgOpt','location','best')





%% Uniform Adaptive System: Mopt, EE and UR 


% figure %8
% hold on; grid on; box on
% plot(MgOpt*ones(1,KgOpt),'-.b')
% plot(Uni_Mopt_Lo(1,:),'.-g')
% plot(Uni_Mopt_Lo(5,:),'.-k')
% plot(Uni_Mopt_Lo(10,:),'.-r')
% xlim([0 110])
% ylim([20 250])
% title('Uniform User Dist')
% legend('Fix MgOpt','Load 10%','Load 50%','Load 100%','location','best')
% ylabel('Number of Active Antennas, M');
% xlabel('Number of active users in the cell, K');

figure %9
hold on; grid on; box on
plot(loading,Uni_AvgEEperBS/1e6,'-xg')
plot(loading,Uni_AvgEEperBS_Fix/1e6,'-xr')
ylim([0 20])
xlim([0 100])
title('Uniform User Dist')
ylabel('Energy Efficiency [Mbits/Joule]');
xlabel('Network load (%)')
legend('Adaptive','Fix MgOpt','location','best')

figure %10
hold on; grid on; box on
plot(loading,Uni_AvgURperBS/1e6,'-xg')
plot(loading,Uni_AvgURperBS_Fix/1e6,'-xr')
ylim([0 250])
xlim([0 100])
title('Uniform User Dist')
ylabel('Average rate per User [Mbps]');
xlabel('Network load (%)');
legend('Adaptive','Fix MgOpt','location','best')




%% CF Adaptive System: Mopt, EE and UR 

% figure %11
% hold on; grid on; box on
% plot(MgOpt*ones(1,KgOpt),'-.b')
% plot(CF_Mopt_Lo(1,:),'.-g')
% plot(CF_Mopt_Lo(5,:),'.-k')
% plot(CF_Mopt_Lo(10,:),'.-r')
% xlim([0 110])
% ylim([20 250])
% title('Center Focused User Dist')
% legend('Fix MgOpt','Load 10%','Load 50%','Load 100%','location','best')
% ylabel('Number of Active Antennas, M');
% xlabel('Number of active users in the cell, K');

figure %12
hold on; grid on; box on
plot(loading,CF_AvgEEperBS/1e6,'-xg')
plot(loading,CF_AvgEEperBS_Fix/1e6,'-xr')
ylim([0 20])
xlim([0 100])
title('Center Focused User Dist')
ylabel('Energy Efficiency [Mbits/Joule]');
xlabel('Network load (%)')
legend('Adaptive','Fix MgOpt','location','best')

figure %13
hold on; grid on; box on
plot(loading,CF_AvgURperBS/1e6,'-xg')
plot(loading,CF_AvgURperBS_Fix/1e6,'-xr')
ylim([0 250])
xlim([0 100])
title('Center Focused User Dist')
ylabel('Average rate per User [Mbps]');
xlabel('Network load (%)');
legend('Adaptive','Fix MgOpt','location','best')











%% EE on BF,Uniform, & CF ULDs 

load('Res_CF.mat')
load('Res_BF.mat')
load('Res_Uniform.mat')


figure  
hold on; box on
plot(loading,CF_AvgEEperBS/1e6,'-^g','linewidth',1.1)
plot(loading,Uni_AvgEEperBS/1e6,'-or','linewidth',1.1)
plot(loading,BF_AvgEEperBS/1e6,'-xb','linewidth',1.1)

plot(loading,CF_AvgEEperBS_Fix/1e6,'--^g','linewidth',1.1)
plot(loading,Uni_AvgEEperBS_Fix/1e6,'--or','linewidth',1.1)
plot(loading,BF_AvgEEperBS_Fix/1e6,'--xb','linewidth',1.1)

ylim([0 20])
xlim([0 100])
legend('Adap. Sys. CF ULD','Adap. Sys. Uni ULD','Adap. Sys. BF ULD','Fixed Sys. CF ULD','Fixed Sys. Uni ULD','Fixed Sys. BF ULD','location','northwest','FontSize',8)
ylabel('Energy Efficiency (Mbits/Joule)');
xlabel('Network load (%)');

set(gcf, 'Position', [500 500 900 500])




%% UR on BF,Uniform, & CF ULDs 
load('Res_CF.mat')
load('Res_BF.mat')
load('Res_Uniform.mat')


figure  
hold on; box on
plot(loading,CF_AvgURperBS_Fix/1e6,'--^g','linewidth',1.1)
plot(loading,CF_AvgURperBS/1e6,'-^g','linewidth',1.1)

plot(loading,Uni_AvgURperBS_Fix/1e6,'--or','linewidth',1.1)
plot(loading,Uni_AvgURperBS/1e6,'-or','linewidth',1.1)

plot(loading,BF_AvgURperBS_Fix/1e6,'--xb','linewidth',1.1)
plot(loading,BF_AvgURperBS/1e6,'-xb','linewidth',1.1)

ylim([0 270])
xlim([0 100])

legend('Fixed Sys. CF ULD','Adap. Sys. CF ULD','Fixed Sys. Uni ULD','Adap. Sys. Uni. ULD','Fixed Sys. BF ULD','Adap. Sys. BF ULD','location','northeast','FontSize',8)

ylabel('Average User-rate (Mbps)');
xlabel('Network load (%)');

set(gcf, 'Position', [500 500 900 500])




























%% Paper Figures

load('UsersDist.mat');

figure
hold on; box on;
% set(gca,'fontsize',12)
% title('User Distribution while serving different cell load')
plot(1:KgOpt,Pi_10,'.-g',1:KgOpt,Pi_50,'.-b',1:KgOpt,Pi_100,'.-r','linewidth',1.3);
xlabel('Number of Active users (K)')
ylabel('User Steady State Probability')
ylim([0 0.25])
legend('10% load','50% load','100% load', 'location', 'best')
set(gcf, 'Position', [500 500 700 400])


%%



load('Res_BF.mat');
load('Res_Uniform.mat');
load('Res_CF.mat');


figure %15
hold on; box on
plot(MgOpt*ones(1,KgOpt),'-.b','linewidth',1.3)

%%%%%%%

plot(CF_Mopt_Lo(1,:),'.-g','linewidth',1.1)
plot(CF_Mopt_Lo(5,:),'.-k','linewidth',1.1)
plot(CF_Mopt_Lo(10,:),'.-r','linewidth',1.1)

annotation(gcf,'ellipse',...
    [0.766471449487556 0.471988795518207 0.0161674394013329 0.042016806722689]);

annotation(gcf,'textarrow',[0.797591055786077 0.772821576763486],...
    [0.449823407623921 0.498964803312629],'TextEdgeColor',[0 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Center-','Focused'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12);

%%%%%%%

plot(BF_Mopt_Lo(1,:),'v-g','linewidth',1.1)
plot(BF_Mopt_Lo(5,:),'v-k','linewidth',1.1)
plot(BF_Mopt_Lo(10,:),'v-r','linewidth',1.1)

annotation(gcf,'ellipse',...
    [0.40625 0.505602240896359 0.0444444444444445 0.301120448179272]);

annotation(gcf,'textarrow',[0.454538062644061 0.429166666666667],...
    [0.678176673895367 0.658263305322129],'TextEdgeColor',[0 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Boundary-','Focused'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12);

%%%%%%%

plot(Uni_Mopt_Lo(1,:),'s-g','linewidth',1.1)
plot(Uni_Mopt_Lo(5,:),'s-k','linewidth',1.1)
plot(Uni_Mopt_Lo(10,:),'s-r','linewidth',1.1)

annotation(gcf,'ellipse',...
    [0.65625 0.547619047619048 0.0319444444444444 0.170868347338936]);

annotation(gcf,'textarrow',[0.698372329184709 0.675],...
    [0.631987944061004 0.620448179271709],'TextEdgeColor',[0 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Uniform'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12);

%%%%%%%

xlim([0 100])
ylim([20 250])
legend('Fix MgOpt','BF, Load 10%','BF, Load 50%','BF, Load 100%',...
       'Uniform, Load 10%','Uniform, Load 50%','Uniform, Load 100%',...
       'CF, Load 10%','CF, Load 50%','CF, Load 100%',...
       'location','best')
 
% legend1 = legend('Fix MgOpt','Cell Load 10%','Cell Load 50%','Cell Load 100%');
% set(legend1,'FontWeight','bold','FontSize',14,'location','best')

   ylabel('Optimal Number of Active Antennas, M','FontWeight','bold','FontSize',14);
xlabel('Number of active users in the cell, K','FontWeight','bold','FontSize',14);
set(gca,'fontsize',12)


%% 


BF_EE_Gain = 100*(BF_AvgEEperBS./BF_AvgEEperBS_Fix-1);
disp(['Average BF EE Gain = ' num2str(mean(BF_EE_Gain))]);
Uni_EE_Gain = 100*(Uni_AvgEEperBS./Uni_AvgEEperBS_Fix-1);
disp(['Average Uni EE Gain = ' num2str(mean(Uni_EE_Gain))]);
CF_EE_Gain = 100*(CF_AvgEEperBS./CF_AvgEEperBS_Fix-1);
disp(['Average CF EE Gain = ' num2str(mean(CF_EE_Gain))]);

figure %16
hold on;box on
plot(loading,BF_EE_Gain,'-.xb','linewidth',1.3)
plot(loading,Uni_EE_Gain,'-.or','linewidth',1.3)
plot(loading,CF_EE_Gain,'-.*g','linewidth',1.3)
xlim([0 100])
ylim([0 300])

ylabel('EE Gain (%)');
xlabel('Network load (%)');
legend('BF ULD','Uni ULD','CF ULD','location','best')

set(gcf, 'Position', [500 500 900 500])



%% 

BF_UR_Loss = 100*(1 - BF_AvgURperBS./BF_AvgURperBS_Fix);
disp(['Average BF UR loss = ' num2str(mean(BF_UR_Loss))]);
Uni_UR_Loss = 100*(1 - Uni_AvgURperBS./Uni_AvgURperBS_Fix);
disp(['Average Uni UR loss = ' num2str(mean(Uni_UR_Loss))]);
CF_UR_Loss = 100*(1 - CF_AvgURperBS./CF_AvgURperBS_Fix);
disp(['Average CF UR loss = ' num2str(mean(CF_UR_Loss))]);

figure %17
hold on;  box on
plot(loading,BF_UR_Loss,'-.xb','linewidth',1.3)
plot(loading,Uni_UR_Loss,'-.or','linewidth',1.3)
plot(loading,CF_UR_Loss,'-.*g','linewidth',1.3)
xlim([0 100])
ylim([0 50])
ylabel('User Rate Loss (%)');
xlabel('Network load (%)');
legend('BF ULD','Uni ULD','CF ULD','location','best')
set(gcf, 'Position', [500 500 900 500])


%% Average number of antennas Vs. Loading for different ULDs

load('Res_BF.mat')
load('Res_CF.mat')
load('Res_Uniform.mat')

figure
hold on; box on
plot(loading,MgOpt*ones(1,10),'-.k','linewidth',1.3)
plot(loading,BF_Mavg_Lo,'-.xr','linewidth',1.3)
plot(loading,Uni_Mavg_Lo,'-.ob','linewidth',1.3)
plot(loading,CF_Mavg_Lo,'-.*g','linewidth',1.3)
xlim([0 100])
ylim([0 250])
ylabel('Average Number of Antennas');
xlabel('Network load (%)');
legend('Fixed Sys. BF ULD','Adap. Sys. BF ULD','Adap. Sys. Uni ULD','Adap. Sys. CF ULD','location','best')

set(gcf, 'Position', [500 500 900 500])


BF_coefficients = polyfit(loading, BF_Mavg_Lo, 1);
BF_slope = BF_coefficients(1)

Uni_coefficients = polyfit(loading, Uni_Mavg_Lo, 1);
Uni_slope = Uni_coefficients(1)

CF_coefficients = polyfit(loading, CF_Mavg_Lo, 1);
CF_slope = CF_coefficients(1)




%% Average Number of Antennas Vs. Loading for different ULDs

% load('Res_BF.mat')
% load('Res_CF.mat')
% load('Res_Uniform.mat')
% 
% figure
% hold on; grid on; box on
% 
% %plot(MgOpt*ones(1,KgOpt),'-.k','linewidth',1.3)
% 
% plot(BF_Mavg_Lo,BF_AvgURperBS/1e6,'-.xr','linewidth',1.3)
% plot(Uni_Mavg_Lo,Uni_AvgURperBS/1e6,'-.ob','linewidth',1.3)
% plot(CF_Mavg_Lo,CF_AvgURperBS/1e6,'-.*g','linewidth',1.3)
% xlim([0 230])
% ylim([0 160])
% xlabel('Average Number of Antennas');
% ylabel('Average User-rate [Mbps]');
% legend('Boundary-focused ULD','Uniform ULD','Center-focused ULD','location','best')



%% Average User-rate Vs. Avrage Number of Antennas

% load('Res_BF.mat')
% load('Res_CF.mat')
% load('Res_Uniform.mat')
% 
% figure
% hold on; grid on; box on
% 
% plot(BF_AvgURperBS/1e6,BF_Mavg_Lo,'-.xr','linewidth',1.3)
% plot(Uni_AvgURperBS/1e6,Uni_Mavg_Lo,'-.ob','linewidth',1.3)
% plot(CF_AvgURperBS/1e6,CF_Mavg_Lo,'-.*g','linewidth',1.3)
% ylim([0 230])
% xlim([0 160])
% ylabel('Average Number of Antennas');
% xlabel('Average User-rate [Mbps]');
% legend('Boundary-focused ULD','Uniform ULD','Center-focused ULD','location','best')



%% Average number of antennas Vs. Loading for BF ULDs

load('Res_BF.mat')
figure
hold on; grid on; box on
plot(MgOpt*ones(1,KgOpt),'-.k','linewidth',1.3)
plot(loading,BF_Mavg_Lo,'-.xr','linewidth',1.3)
xlim([0 100])
ylim([0 250])
ylabel('Average Number of Antennas');
xlabel('Network load (%)');
legend('BF, M_{Fix}','BF, M_{Opt}','location','best')


%% UR on BF ULDs 
load('Res_BF.mat');
figure  
hold on; grid on; box on
plot(loading,BF_AvgURperBS_Fix/1e6,'-.k','linewidth',1.1)
plot(loading,BF_AvgURperBS/1e6,'-.xr','linewidth',1.1)
ylim([0 270])
xlim([0 100])
legend('BF, M_{Fix}','BF, M_{Opt}','location','best')
ylabel('Average User-rate [Mbps]');
xlabel('Network load (%)');





















%% UR on Uniform ULDs
load('Res_Uniform.mat');
figure  
hold on; grid on; box on
plot(loading,Uni_AvgURperBS_Fix/1e6,'--or','linewidth',1.1)
plot(loading,Uni_AvgURperBS/1e6,'--ob','linewidth',1.1)
ylim([0 270])
xlim([0 100])
legend('Uniform, M_{Fix}','Uniform, M_{Opt}','location','best')
ylabel('Average User-rate [Mbps]');
xlabel('Network load (%)');
%% UR on CF ULDs
load('Res_CF.mat');
figure
hold on; grid on; box on
plot(loading,CF_AvgURperBS_Fix/1e6,'-xr','linewidth',1.1)
plot(loading,CF_AvgURperBS/1e6,'-xb','linewidth',1.1)
ylim([0 270])
xlim([0 100])
legend('CF, M_{Fix}','CF, M_{Opt}','location','best')
ylabel('Average User-rate [Mbps]');
xlabel('Network load (%)');

