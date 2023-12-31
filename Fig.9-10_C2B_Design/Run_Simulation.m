%% Simulation Enviroment
clear all; close all; clc
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

p_mean_per_AE = 20/222;

Mmin = ceil(Pc/p_max);
GoS = 0.02;

%%
% PLO PLI Calculation 
Ri = [Rmin 400 Rmax];
Cen2Bou_Ratio = [1:-.1:.8 0.7725 0.7:-.1:0; 0:.1:.2 0.2275 0.3:.1:1];
PLO_Network_BSs = cell(12,1);
PLI_Network_BSs = cell(12,1);
UELocations_C2B = cell(1,1);
for r = 1:12
    UELocations_C2B{r} = UE_insertion_MonteCarlo_HexCell(TestPoints,Cen2Bou_Ratio(:,r)',Ri,false);
    [PLO_Network_BSs{r}, PLI_Network_BSs{r}] = Wrap_Around_PLO_PLI(BSLocations,UELocations_C2B{r},1,Rmax,false);
end


%%

MgOpt = zeros(1,12);
KgOpt = zeros(1,12);
K_design = zeros(1,12);
M_design = zeros(1,12);
lambdaS = zeros(1,12);
Mopt = cell(1,12);
AvgEEperBS = zeros(1,12);
AvgURperBS = zeros(1,12);

Mmax = 300;Kmax=150;
for r = 1:12
    % Fixed System Dimensioning at 100% cell loading
    disp(['Current Cen2Bou Ratio = ' num2str(Cen2Bou_Ratio(:,r)') ' %'])    
    [MgOpt(r),KgOpt(r)] = Global_Local_Opt_EE(PLO_Network_BSs{r},PLI_Network_BSs{r},Kmax,Mmax,Pc,PmaxPA,false);
    disp(['KgOpt = ' num2str(KgOpt(r)) '  MgOpt = ' num2str(MgOpt(r))])
        
%     % Fixed System Dimensioning at 100% cell loading (NO NEED FOR THIS STEP)
%     Rc = zeros(1,KgOpt(r));
%     Mc = MgOpt(r); Md = MgOpt(r)*ones(1,18);p=Pc/MgOpt(r);
%     for K= 1:KgOpt(r)
%         [~,Rc(K),~] = EE_R_Ptot_PA(PLO_Network_BSs{r,1},PLI_Network_BSs{r,1},KgOpt(r),K,Mc,Md,p,PmaxPA,[]);
%     end
%     lambdaS(r) = Searching_lambdaS(Rc,KgOpt(r),GoS,false);
%     [Mopt_temp,AvgEEperBS(r), AvgURperBS(r)] = SysOpt_Adaptive_M_SingleCell(PLO_Network_BSs{r},PLI_Network_BSs{r},KgOpt(r),Mmin,MgOpt(r),Pc,PmaxPA,100,lambdaS(r));
%     K_design(r) = length(Mopt_temp);
%     M_design(r) = Mopt_temp(end);
%     Mopt{r}=Mopt_temp;
%     disp(['K_design = ' num2str(K_design(r)) '  M_design = ' num2str(M_design(r))])
    
end

%%

EE_MK_design = zeros(1,12);
Rt = zeros(1,12);
Rc = zeros(1,12);
Rb = zeros(1,12);

for r =1:12
    Mc = M_design(r);
    Md = Mc*ones(1,18);
    p=Pc/Mc;
   [EE_MK_design(r),~,~,~,R_TP,~] = EE_R_Ptot_PA(PLO_Network_BSs{r},PLI_Network_BSs{r},K_design(r),K_design(r),Mc,Md,p,PmaxPA,[]);
   Rt(r) = mean(R_TP);
   UE_b = abs(UELocations_C2B{r})>=400;
   UE_c = ~ UE_b;
   Rc(r) = mean(R_TP(UE_c));
   Rb(r) = mean(R_TP(UE_b));
end

%%
save('Res_C2B_Design.mat','Cen2Bou_Ratio','KgOpt','MgOpt','Rb','Rc','Rt','EE_MK_design');





