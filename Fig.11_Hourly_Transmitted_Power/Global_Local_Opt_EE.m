function [MgOpt,KgOpt,Mopt,EEopt,Ropt] = Global_Local_Opt_EE(PLO,PLI,Kmax,Mmax,Pc,PmaxPA,plott)

%%
% clear all; close all; clc
% % Simulation Enviroment
% Rmax = 500; %Cell radius (distance to a vertex of hexagonal cell)
% Rmin = 35; %Users inside this circle will not considered in simulations
% TestPoints = 10000;
% ISD = Rmax*sqrt(3);
% % % Coordinates of all BSs in the Hexagonal Network
% u = [0 1 0 -1 -1  0  1  2 2 1 0 -1 -2 -2 -2 -1  0  1  2]; % 30-Degree axis
% v = [0 0 1  1  0 -1 -1 -1 0 1 2  2  2  1  0 -1 -2 -2 -2]; % Vertical axis
% % Setting the BSs in their locations as seen from the origin.
% BSLocations = sqrt(3).*(ISD/2+1i*Rmax/2).*u + (0+1i*ISD).*v;
% Pc = 20;
% Mmax = 300;Kmax=150;
% PmaxPAcriteria = 'Fixed';
% PmaxPA = 10^0.6; % 6dB
% GoS = 0.02;
% % EE Vs. Boundary User Rate
% UELocations = cell(1,9);
% PLO_Network = cell(1,19);
% PLI_Network = cell(1,19);
% Ri = [Rmin 400 Rmax];
% Cen2Bou_Ratio = [1:-.1:0; 0:.1:1];
% for w = 1:11
% UELocations{w} = UE_insertion_MonteCarlo_HexCell(TestPoints,Cen2Bou_Ratio(:,w),Ri,false);
% [PLO_Network{w},PLI_Network{w}] = Wrap_Around_PLO_PLI(BSLocations,UELocations{w},1,Rmax,false);
% end
% PLO = PLO_Network{4};
% PLI = PLI_Network{4};
% plott  = true;
% 


% This function search for the optimal Operating Point when BSs are
% transmitting constant average power, Pc






EE = zeros(Mmax,Kmax);
R = zeros(Mmax,Kmax);
p_max = PmaxPA/10^0.8;
Mmin = ceil(Pc/p_max);
for K=1:Kmax
    %disp(['Current K = ' num2str(K) '   Kmax = ' num2str(Kmax)]);
    for M=Mmin:Mmax;
        p = Pc/M;
         [EE(M,K),R(M,K),~,~] = EE_R_Ptot_PA(PLO,PLI,K,K,M,M*ones(1,18),p,PmaxPA,[]);
    end
end

[EEmax,Mind]  = max(EE);
[EEgOpt, KgOpt] = max(EEmax);
MgOpt = Mind(KgOpt);

Mopt = Mind(1:KgOpt);
EEopt = EEmax(1:KgOpt);

Rind = zeros(1,Kmax);
for K=1:Kmax
    M = Mind(K);
    Rind(K) = R(M,K);
end
Ropt = Rind(1:KgOpt);




if plott ==true
    
    % Illustrating Global Optimum point at EE
    % Energy Efficiency [Mbit/Joule]
    figure
    hold on; box off; grid on;
    title('Energy Efficiency [Mbit/Joule]')
    gridDensity = 25;
    surface(1:Kmax,1:Mmax,EE/1e6,'EdgeColor','none'); %Plot the 3d surface
    colormap(autumn);
    
    %Plot lines on top of the 3d surface, to make it easier to see the shape
    for m = [1 gridDensity:gridDensity:Mmax]
        plot3(1:Kmax,m*ones(1,Kmax),EE(m,:)/1e6,'k-');
    end
    for k = [1 gridDensity:gridDensity:Kmax]
        plot3(k*ones(1,Mmax),1:Mmax,EE(:,k)/1e6,'k-');
    end
    
    %Compute and plot the optimal point 
    plot3(KgOpt,MgOpt,EEgOpt/1e6,'b*','MarkerSize',13);    
    txt = {'Global Optimum:',['M = ' num2str(MgOpt) ', K = ' num2str(KgOpt)],['EE = ' num2str(EEgOpt/1e6) ' Mbit/J']};
    text(Kmax,100,EEgOpt/1e6,txt,'HorizontalAlignment','left','BackgroundColor', [1 1 1],'EdgeColor','b')
    view([-46 24]);
    axis([0 Kmax 0 Mmax]);
    xlabel('Number of Users (K)');
    ylabel('Number of Antennas (M)');
    zlabel('Energy Efficiency [Mbit/Joule]');
    
    % % Optimal Operating Point for 1:KgOpt
    figure
    subplot(3,1,1);
    hold on; grid on
    plot(1:Kmax,Mind)
    plot(KgOpt,MgOpt,'ro')
    ylabel('Optimal M')
    
    subplot(3,1,2)
    hold on; grid on
    plot(1:Kmax,EEmax/1e6)
    plot(KgOpt,EEgOpt/1e6,'ro')
    ylabel('Max EE')

    subplot(3,1,3);
    hold on; grid on
    plot(1:Kmax,Rind/1e6)
    plot(KgOpt,Rind(KgOpt)/1e6,'ro')
    xlabel('Number of Users,K')
    ylabel('User Rate')

    
end

end