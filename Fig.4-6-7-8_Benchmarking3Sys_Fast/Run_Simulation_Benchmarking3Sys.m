clear all,close all;clc

%% Simulation Enviroment
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






%% Dimensioning Fixed Antenna System at BF ULD
disp('Searching for Global Optimum E at BF ULD');
Ri = [Rmin 200 400 Rmax];
ULD_BF = [0.10 0.20 0.70]';
UELocations = UE_insertion_MonteCarlo_HexCell(TestPoints,ULD_BF,Ri,true);


Ni = 1;
PLO_Network = cell(1,Ni);
PLI_Network = cell(1,Ni);
for BSj=1:Ni
    [PLO_Network{BSj},  PLI_Network{BSj}] = Wrap_Around_PLO_PLI(BSLocations,UELocations,BSj,Rmax,true);
end


%[MgOpt,KgOpt] = Global_Local_Opt_EE(PLO_Network{1},PLI_Network{1},Kmax,Mmax,Pc,PmaxPA,false);
KgOpt = 96;  MgOpt = 228; 
disp(['KgOpt = ' num2str(KgOpt) ' ,  MgOpt = ' num2str(MgOpt)])


% % Dimensioning Fixed Antenna System when Loading 100%
% Rc = zeros(1,KgOpt);
% Mc = MgOpt; Md = MgOpt*ones(1,18);p=Pc/MgOpt;
% for K= 1:KgOpt
%     [~,Rc(K),~] = EE_R_Ptot_PA(PLO_Network{1},PLI_Network{1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
% end
% lambdaS = Searching_lambdaS(Rc,KgOpt,GoS,false);
% [Mopt_BF_CL100] = SysOpt_Adaptive_M_SingleCell(PLO_Network{1},PLI_Network{1},KgOpt,Mmin,MgOpt,Pc,PmaxPA,100,lambdaS);
% MgOpt = Mopt_BF_CL100(end);
% MgOpt = 228;
% disp(['KgOpt = ' num2str(KgOpt) ' ,  Mopt(end) = ' num2str(MgOpt)])

Rc = zeros(1,KgOpt);
Mc = MgOpt; Md = MgOpt*ones(1,18);p=Pc/MgOpt;
for K= 1:KgOpt
    [~,Rc(K),~] = EE_R_Ptot_PA(PLO_Network{1},PLI_Network{1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
end
lambdaS = Searching_lambdaS(Rc,KgOpt,GoS,false);








%% Users Dist at 10,50,100 Cell Loading
% Pi_10 = zeros(1,KgOpt);
% Pi_50 = zeros(1,KgOpt);
% Pi_100 = zeros(1,KgOpt);
% for K=1:KgOpt
%     Pi_10(K) = MGmm_SD_Queue(K,KgOpt,0.10*lambdaS,1,Rc);
%     Pi_50(K) = MGmm_SD_Queue(K,KgOpt,0.50*lambdaS,1,Rc);
%     Pi_100(K) = MGmm_SD_Queue(K,KgOpt,lambdaS,1,Rc);
% end
% save('UsersDist.mat','KgOpt','Pi_10','Pi_50','Pi_100');








%% Performance of Adaptive System in BF, Uniform and CF ULD:
%
disp('$$$$$$$$$$$$$$$$$$$$$ Adaptive System in BF ULD $$$$$$$$$$$$$$$$$$$$$$$$$');
BF_Mopt_Lo = zeros(length(loading),KgOpt);
BF_Mavg_Lo = zeros(1,length(loading));
BF_AvgEEperBS = zeros(1,length(loading));
BF_AvgURperBS = zeros(1,length(loading));
for L = 1:length(loading);
    [BF_Mopt_Lo(L,:),BF_AvgEEperBS(L), BF_AvgURperBS(L),~,BF_Mavg_Lo(L)] = SysOpt_Adaptive_M_SingleCell(PLO_Network{1},PLI_Network{1},KgOpt,Mmin,MgOpt,Pc,PmaxPA,loading(L),lambdaS);
end

BF_AvgEEperBS_Fix = zeros(1,length(loading));
BF_AvgURperBS_Fix = zeros(1,length(loading));
for L = 1:length(loading);
    [BF_AvgEEperBS_Fix(1,L), BF_AvgURperBS_Fix(1,L)] = SysOpt_Fixed_M(PLO_Network,PLI_Network,KgOpt,MgOpt,Pc,PmaxPA,loading(L),lambdaS);
end

save('Res_BF.mat','BF_Mopt_Lo','BF_Mavg_Lo','BF_AvgEEperBS','BF_AvgURperBS','BF_AvgEEperBS_Fix','BF_AvgURperBS_Fix','KgOpt','MgOpt')










%%
disp('$$$$$$$$$$$$$$$$$$ Adaptive System in Uniform ULD $$$$$$$$$$$$$$$$$$$$$$$');
Ri = [Rmin Rmax];
ULD_Uni = 1 ;
UELocations = UE_insertion_MonteCarlo_HexCell(TestPoints,ULD_Uni,Ri,false);
PLO_Network = cell(1,19);
PLI_Network = cell(1,19);
for BSj=1:19
[PLO_Network{BSj},  PLI_Network{BSj}] = Wrap_Around_PLO_PLI(BSLocations,UELocations,BSj,Rmax,false);
end

Uni_Mopt_Lo = zeros(length(loading),KgOpt);
Uni_Mavg_Lo = zeros(1,length(loading));
Uni_AvgEEperBS = zeros(1,length(loading));
Uni_AvgURperBS = zeros(1,length(loading));
for L = 1:length(loading);
    [Uni_Mopt_Lo(L,:),Uni_AvgEEperBS(L), Uni_AvgURperBS(L),~,Uni_Mavg_Lo(L)] = SysOpt_Adaptive_M_SingleCell(PLO_Network{1},PLI_Network{1},KgOpt,Mmin,MgOpt,Pc,PmaxPA,loading(L),lambdaS);
end
Uni_AvgEEperBS_Fix = zeros(1,length(loading));
Uni_AvgURperBS_Fix = zeros(1,length(loading));
for L = 1:length(loading);
    [Uni_AvgEEperBS_Fix(1,L), Uni_AvgURperBS_Fix(1,L)] = SysOpt_Fixed_M(PLO_Network,PLI_Network,KgOpt,MgOpt,Pc,PmaxPA,loading(L),lambdaS);
end
save('Res_Uniform.mat','Uni_Mopt_Lo','Uni_Mavg_Lo','Uni_AvgEEperBS','Uni_AvgURperBS','Uni_AvgEEperBS_Fix','Uni_AvgURperBS_Fix','KgOpt','MgOpt');













%%
disp('$$$$$$$$$$$$$$$$$$$$$ Adaptive System in CF ULD $$$$$$$$$$$$$$$$$$$$$$$$$');
Ri = [Rmin 200 400 Rmax];
ULD_CF = [0.80 0.10 0.10]';
UELocations = UE_insertion_MonteCarlo_HexCell(TestPoints,ULD_CF,Ri,true);
PLO_Network = cell(1,19);
PLI_Network = cell(1,19);
for BSj=1:19
[PLO_Network{BSj},PLI_Network{BSj}] = Wrap_Around_PLO_PLI(BSLocations,UELocations,BSj,Rmax,false);
end

CF_Mopt_Lo = zeros(length(loading),KgOpt);
CF_Mavg_Lo = zeros(1,length(loading));
CF_AvgEEperBS = zeros(1,length(loading));
CF_AvgURperBS = zeros(1,length(loading));
for L = 1:length(loading);
    [CF_Mopt_Lo(L,:),CF_AvgEEperBS(L), CF_AvgURperBS(L),~,CF_Mavg_Lo(L)] = SysOpt_Adaptive_M_SingleCell(PLO_Network{1},PLI_Network{1},KgOpt,Mmin,MgOpt,Pc,PmaxPA,loading(L),lambdaS);
end

CF_AvgEEperBS_Fix = zeros(1,length(loading));
CF_AvgURperBS_Fix = zeros(1,length(loading));
for L = 1:length(loading);
    [CF_AvgEEperBS_Fix(1,L), CF_AvgURperBS_Fix(1,L)] = SysOpt_Fixed_M(PLO_Network,PLI_Network,KgOpt,MgOpt,Pc,PmaxPA,loading(L),lambdaS);
end
save('Res_CF.mat','CF_Mopt_Lo','CF_Mavg_Lo','CF_AvgEEperBS','CF_AvgURperBS','CF_AvgEEperBS_Fix','CF_AvgURperBS_Fix','KgOpt','MgOpt');





