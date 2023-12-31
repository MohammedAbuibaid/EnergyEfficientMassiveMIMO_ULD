%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scenario A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Enviroment
clear all,close all;clc
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
PmaxPA = 10^0.6; % 6dB
p_max = PmaxPA/10^0.8;
Mmin = ceil(Pc/p_max);
GoS = 0.02;
Ri = [Rmin 200 400 Rmax];


%% ULD Model
disp('%%%%%%%%%%%%%%%%%%%%%%% Scenario A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%    Time:Hours    01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21  22  23 24
Loading = [        63 44 25 19 13 13 19 25 38 50 56 63 63 69 75 75 75 81 88 94 100 100 94 81];
C2M2B_Scenario1 = [10 10 10 10 10 10 10 60 60 80 80 80 80 80 80 80 80 80 80 60 60  10  10 10;
                   20 20 20 20 20 20 20 20 20 10 10 10 10 10 10 10 10 10 10 20 20  20  20 20;
                   70 70 70 70 70 70 70 20 20 10 10 10 10 10 10 10 10 10 10 20 20  70  70 70]./100;                   
PLO_Network_hr = cell(24,19);
PLI_Network_hr = cell(24,19);
UELocations_hr = cell(24,1);
for h = 1:24
    UELocations_hr{h} = UE_insertion_MonteCarlo_HexCell(TestPoints,C2M2B_Scenario1(:,h)',Ri,false);
    for BSj=1:19
        [PLO_Network_hr{h,BSj}, PLI_Network_hr{h,BSj}] = Wrap_Around_PLO_PLI(BSLocations,UELocations_hr{h},BSj,Rmax,false);
    end
end


%%
% Dimensioning Fixed Antenna System at BF hour h = 22
disp('Dimensioning Fixed Antenna System at hour h = 22');
h_design = 22;
disp(['Current Cen2Bou Ratio = ' num2str(C2M2B_Scenario1(:,h_design)') ' %'])
Mmax = 250;Kmax=125;
[MgOpt,KgOpt] = Global_Local_Opt_EE(PLO_Network_hr{h_design,1},PLI_Network_hr{h_design,1},Kmax,Mmax,Pc,PmaxPA,false);
% KgOpt = 97; MgOpt = 228;
disp(['KgOpt = ' num2str(KgOpt)   ' MgOpt= ' num2str(MgOpt)])


% % Dimensioning Fixed Antenna System when Loading 100%
% Rc = zeros(1,KgOpt);
% Mc = MgOpt; Md = MgOpt*ones(1,18); p=Pc/MgOpt;
% for K= 1:KgOpt
%     [~,Rc(K),~] = EE_R_Ptot_PA(PLO_Network_hr{h_design,1},PLI_Network_hr{h_design,1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
% end
% lambdaS = Searching_lambdaS(Rc,KgOpt,GoS,false);
% [Mopt_h22] = SysOpt_Adaptive_M_SingleCell(PLO_Network_hr{h_design,1},PLI_Network_hr{h_design,1},KgOpt,Mmin,MgOpt,Pc,PmaxPA,100,lambdaS);
% K_design = length(Mopt_h22); M_design = Mopt_h22(end);
% % K_design = 97 ;  M_design = 228;
% disp(['K_design = ' num2str(K_design) ' ,  M_design = ' num2str(M_design)])
% % Update KgOpt and MgOpt
% MgOpt = M_design;
% KgOpt = K_design;

Rc = zeros(1,KgOpt);
Mc = MgOpt; Md = MgOpt*ones(1,18);p=Pc/MgOpt;
for K= 1:KgOpt
    [~,Rc(K),~] = EE_R_Ptot_PA(PLO_Network_hr{h_design,1},PLI_Network_hr{h_design,1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
end
lambdaS = Searching_lambdaS(Rc,KgOpt,GoS,false);


%% Performance of Adaptive System according to Scenario A:
Mopt= cell(24,1);
AvgEEperBS= zeros(24,1);
AvgURperBS = zeros(24,1);
AvgPtotperBS = zeros(24,1);
for h = 1:24
    disp(['****** Adaptive Sys ****** Scenario A ****** Time ( Hour = ' num2str(h)  ' ) ******']);
    disp(['Loading = ' num2str(Loading(h))])
    disp(['Current C2M2B Ratio = ' num2str(C2M2B_Scenario1(:,h)') ' %'])
    [Mopt{h},AvgEEperBS(h),AvgURperBS(h),AvgPtotperBS(h)] = SysOpt_Adaptive_M_SingleCell(PLO_Network_hr{h,1},PLI_Network_hr{h,1},KgOpt,Mmin,MgOpt,Pc,PmaxPA,Loading(h),lambdaS);
end

AvgEEperBS_Fix = zeros(24,1);
AvgURperBS_Fix = zeros(24,1);
AvgPtotperBS_Fix = zeros(24,1);
for h = 1:24
    disp(['****** Fixed Sys ****** Scenario A ****** Time ( Hour = ' num2str(h)  ' ) ******']);
    disp(['Loading = ' num2str(Loading(h))])
    disp(['Current C2M2B Ratio = ' num2str(C2M2B_Scenario1(:,h)') ' %'])
    [AvgEEperBS_Fix(h), AvgURperBS_Fix(h),AvgPtotperBS_Fix(h)] = SysOpt_Fixed_M(PLO_Network_hr(h,:),PLI_Network_hr(h,:),KgOpt,MgOpt,Pc,PmaxPA,Loading(h),lambdaS);
end

%%

Mmax_h = zeros(1,24);
WA_Mopt_h = zeros(1,24);
for h = 1:24
    Mopt_h  = Mopt{h};
    Mmax_h(h) = Mopt_h(end);
    Ropt = zeros(1,KgOpt);
    for K=1:KgOpt
       Mc = Mopt_h(K);
       Md = Mc*ones(1,18);
       p = Pc/Mc;
       [~,Ropt(K),~] = EE_R_Ptot_PA(PLO_Network_hr{h,1},PLI_Network_hr{h,1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
    end
    Pi = zeros(1,KgOpt);
    for K=1:KgOpt
        [Pi(K),~] = MGmm_SD_Queue(K,KgOpt,(Loading(h)/100)*lambdaS,1,Ropt);
    end
    WA_Mopt_h(h) = sum(Mopt_h.*Pi);
end
save('Res_Scenario_A.mat','KgOpt','MgOpt','Mopt','AvgEEperBS','AvgURperBS','AvgPtotperBS','AvgEEperBS_Fix','AvgURperBS_Fix','AvgPtotperBS_Fix','Mmax_h','WA_Mopt_h');
%%


















%% %%%%%%%%%%%%%%%%%%%%%%%%%% Scenario B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%% Scenario B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% ULD Model
%    Time:Hours    01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21  22  23 24
Loading = [        63 44 25 19 13 13 19 25 38 50 56 63 63 69 75 75 75 81 88 94 100 100 94 81];
C2M2B_Scenario2 = [70 70 70 70 70 70 70 20 20 10 10 10 10 10 10 10 10 10 10 20 20  70  70 70;
                   20 20 20 20 20 20 20 20 20 10 10 10 10 10 10 10 10 10 10 20 20  20  20 20;
                   10 10 10 10 10 10 10 60 60 80 80 80 80 80 80 80 80 80 80 60 60  10  10 10]./100;             
PLO_Network_hr = cell(24,19);
PLI_Network_hr = cell(24,19);
UELocations_hr = cell(24,1);
for h = 1:24
    UELocations_hr{h} = UE_insertion_MonteCarlo_HexCell(TestPoints,C2M2B_Scenario2(:,h)',Ri,false);
    for BSj=1:19
        [PLO_Network_hr{h,BSj}, PLI_Network_hr{h,BSj}] = Wrap_Around_PLO_PLI(BSLocations,UELocations_hr{h},BSj,Rmax,false);
    end
end

% Dimensioning Fixed Antenna System at BF hour h = 19
h_design = 19;
disp('Dimensioning Fixed Antenna System at hour h = 19');
disp(['Current Cen2Bou Ratio = ' num2str(C2M2B_Scenario2(:,h_design)') ' %'])
Mmax = 250;Kmax=125;
[MgOpt,KgOpt] = Global_Local_Opt_EE(PLO_Network_hr{h_design,1},PLI_Network_hr{h_design,1},Kmax,Mmax,Pc,PmaxPA,false);
% KgOpt = 111; MgOpt = 167;
disp(['KgOpt = ' num2str(KgOpt)   ' MgOpt = ' num2str(MgOpt)])


% % Dimensioning Fixed Antenna System when Loading 100%
% Rc = zeros(1,KgOpt);
% Mc = MgOpt; Md = MgOpt*ones(1,18);p=Pc/MgOpt;
% for K= 1:KgOpt
%     [~,Rc(K),~] = EE_R_Ptot_PA(PLO_Network_hr{h_design,1},PLI_Network_hr{h_design,1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
% end
% lambdaS = Searching_lambdaS(Rc,KgOpt,GoS,false);
% [Mopt_h19] = SysOpt_Adaptive_M_SingleCell(PLO_Network_hr{h_design,1},PLI_Network_hr{h_design,1},KgOpt,Mmin,MgOpt,Pc,PmaxPA,100,lambdaS);
% K_design = length(Mopt_h19); M_design = Mopt_h19(end);
% disp(['K_design = ' num2str(K_design) ' ,  M_design = ' num2str(M_design)])
% % Update KgOpt and MgOpt
% MgOpt = M_design;
% KgOpt = K_design;

Rc = zeros(1,KgOpt);
Mc = MgOpt; Md = MgOpt*ones(1,18);p=Pc/MgOpt;
for K= 1:KgOpt
    [~,Rc(K),~] = EE_R_Ptot_PA(PLO_Network_hr{h_design,1},PLI_Network_hr{h_design,1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
end
lambdaS = Searching_lambdaS(Rc,KgOpt,GoS,false);


%% Performance of Adaptive System according to Scenario B:
Mopt= cell(24,1);
AvgEEperBS= zeros(24,1);
AvgURperBS = zeros(24,1);
AvgPtotperBS = zeros(24,1);

for h = 1:24
    disp(['****** Adaptive Sys ****** Scenario B ****** Time ( Hour = ' num2str(h)  ' ) ******']);
    disp(['Loading = ' num2str(Loading(h))])
    disp(['Current C2M2B Ratio = ' num2str(C2M2B_Scenario2(:,h)') ' %'])
    [Mopt{h},AvgEEperBS(h),AvgURperBS(h),AvgPtotperBS(h)] = SysOpt_Adaptive_M_SingleCell(PLO_Network_hr{h,1},PLI_Network_hr{h,1},KgOpt,Mmin,MgOpt,Pc,PmaxPA,Loading(h),lambdaS);
end

AvgEEperBS_Fix = zeros(24,1);
AvgURperBS_Fix = zeros(24,1);
AvgPtotperBS_Fix = zeros(24,1);
for h = 1:24
    disp(['****** Fixed Sys ****** Scenario B ****** Time ( Hour = ' num2str(h)  ' ) ******']);
    disp(['Loading = ' num2str(Loading(h))])
    disp(['Current C2M2B Ratio = ' num2str(C2M2B_Scenario2(:,h)') ' %'])
    [AvgEEperBS_Fix(h), AvgURperBS_Fix(h),AvgPtotperBS_Fix(h)] = SysOpt_Fixed_M(PLO_Network_hr(h,:),PLI_Network_hr(h,:),KgOpt,MgOpt,Pc,PmaxPA,Loading(h),lambdaS);
end

%%
Mmax_h = zeros(1,24);
WA_Mopt_h = zeros(1,24);
for h = 1:24
    Mopt_h  = Mopt{h};
    Mmax_h(h) = Mopt_h(end);
    Ropt = zeros(1,KgOpt);
    for K=1:KgOpt
       Mc = Mopt_h(K);
       Md = Mc*ones(1,18);
       p = Pc/Mc;
       [~,Ropt(K),~] = EE_R_Ptot_PA(PLO_Network_hr{h,1},PLI_Network_hr{h,1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
    end
    Pi = zeros(1,KgOpt);
    for K=1:KgOpt
        [Pi(K),~] = MGmm_SD_Queue(K,KgOpt,(Loading(h)/100)*lambdaS,1,Ropt);
    end
    WA_Mopt_h(h) = sum(Mopt_h.*Pi);
end
save('Res_Scenario_B.mat','KgOpt','MgOpt','Mopt','AvgEEperBS','AvgURperBS','AvgPtotperBS','AvgEEperBS_Fix','AvgURperBS_Fix','AvgPtotperBS_Fix','Mmax_h','WA_Mopt_h');
