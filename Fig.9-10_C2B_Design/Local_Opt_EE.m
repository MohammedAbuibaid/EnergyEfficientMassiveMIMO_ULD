function [EEmax,Mopt,Ropt] = Local_Opt_EE(PLO,PLI,KgOpt,MgOpt,Pc,PmaxPA,plott)
% % This function find the optimal Operating points in the ranges 1:Kmax and 1:Mmax
% % when BSs are transmitting constant average power, Pc
p_max = PmaxPA/10^0.8;
Mmin = ceil(Pc/p_max);

EE = zeros(MgOpt,KgOpt);
R = zeros(MgOpt,KgOpt);
for K=1:KgOpt
    %disp(['Current K = ' num2str(K) '   Kmax = ' num2str(Kmax)]);
    for M=Mmin:MgOpt;
        p = Pc/M;
         [EE(M,K),R(M,K),~,~] = EE_R_Ptot_PA(PLO,PLI,KgOpt,K,M,M*ones(1,18),p,PmaxPA,[]);
    end
end
[EEmax,Mopt]  = max(EE);

Ropt = zeros(1,KgOpt);
for K=1:KgOpt
    M = Mopt(K);
    Ropt(K) = R(M,K);
end

%%
if plott ==true
    % % Optimal Operating Point for 1:KgOpt
    figure
    subplot(3,1,1)
    plot(1:KgOpt,EEmax/1e6)
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