function [EEopt,Ropt,Mopt,KgOpt] = PA_Dimen_gOptEE(PLO,PLI,Kmax,Mmax,Pc,PmaxPAcriteria,Pmax_PA,fxy,plott)
% This function search for the optimal Operating Point when BSs are
% transmitting constant average power, Pc

EE = zeros(Mmax,Kmax);
R = zeros(Mmax,Kmax);
p_max = Pmax_PA/10^0.8;
Mmin = ceil(Pc/p_max);
for K=1:Kmax
    %disp(['Current K = ' num2str(K) '   Kmax = ' num2str(Kmax)]);
    for M=Mmin:Mmax;
        p = Pc/M;
%   function [EE,Rc,Ptot,QoS] = EE_R_Ptot_PA(PLO,PLI,Kmax,K,Mc,Md,p,Pmax_PA,mimR,fxy)
         [EE(M,K),R(M,K),~,~] = EE_R_Ptot_PA(PLO,PLI,K,K,M,M*ones(1,18),p,Pmax_PA,[],fxy);
    end
end

[EEvalues,MInd]  = max(EE);
[EEgOpt, KgOpt] = max(EEvalues);
MgOpt = MInd(KgOpt);
Mopt = MInd(1:KgOpt);

EEopt = zeros(1,KgOpt);
Ropt = zeros(1,KgOpt);
for K=1:KgOpt
    M = Mopt(K);
    EEopt(K) = EE(M,K);
    Ropt(K) = R(M,K);
end

%%
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
    subplot(3,1,1)
    plot(1:KgOpt,EEopt/1e6)
    title(['Max EE can be obtained for K = 1:' num2str(KgOpt) ', [Mbit/Joule]'])
    grid on
    
    subplot(3,1,2);
    plot(1:KgOpt,Mopt)
    title(['Optimal M for K = 1:' num2str(KgOpt)])
    grid on
    
    subplot(3,1,3);
    plot(1:KgOpt,Ropt/1e6)
    title('User Rate for for Optimal M, [Mbps]')
    grid on
    ylim([0 250])
    xlabel('Number of Users,K')
end

end