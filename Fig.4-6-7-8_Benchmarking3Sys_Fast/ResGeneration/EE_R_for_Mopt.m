function [EE,R] = EE_R_for_Mopt(PLO,PLI,KgOpt,Mopt,Pc,PmaxPA,plott)

EE = zeros(1,KgOpt);
R = zeros(1,KgOpt);
for K=1:KgOpt
    %disp(['Current K = ' num2str(K) '   Kmax = ' num2str(Kmax)]);
    M = Mopt(K);
    p = Pc/M;
    [EE(K),R(K),~,~] = EE_R_Ptot_PA(PLO,PLI,KgOpt,K,M,M*ones(1,18),p,PmaxPA,[]);
end


%%
if plott ==true
    % % Optimal Operating Point for 1:KgOpt
    figure
    subplot(3,1,1)
    plot(1:KgOpt,EE/1e6)
    title(['EE can be obtained for K = 1:' num2str(KgOpt) ', [Mbit/Joule]'])
    grid on
    
    subplot(3,1,2);
    plot(1:KgOpt,Mopt)
    title(['Optimal M for K = 1:' num2str(KgOpt)])
    grid on
    
    subplot(3,1,3);
    plot(1:KgOpt,R/1e6)
    title('User Rate for for Optimal M, [Mbps]')
    grid on
    ylim([0 250])
    xlabel('Number of Users,K')
end

end