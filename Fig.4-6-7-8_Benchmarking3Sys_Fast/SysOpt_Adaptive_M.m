function [Mopt,AvgEEperBS, AvgURperBS]=SysOpt_Adaptive_M(PLO_Network,PLI_Network,KgOpt,Mmin,MgOpt,Pc,PmaxPA,loading,lambdaS)

Mref=MgOpt*ones(1,19);
disp(['Mref = ' num2str(Mref) ])
Converged=false;
M_iterate=Mref;
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp('For BSj = 1');
while(~Converged)
    Md_temp= M_iterate;
    Md_temp(1)= [];
    % Finding the optimum antennas for each user state
    EE = zeros(MgOpt,KgOpt);
    UserRate = zeros(MgOpt,KgOpt);
    for K=1:KgOpt
        for M=Mmin:MgOpt;
            Mc = M; Md = Md_temp;
            p = Pc/M;
            [EE(M,K),UserRate(M,K)] = EE_R_Ptot_PA(PLO_Network{1},PLI_Network{1},KgOpt,K,Mc,Md,p,PmaxPA,[]);
        end
    end
    [EEopt, Mopt]=max(EE);
    Ropt = zeros(1,KgOpt);
    for K=1:KgOpt
        Ropt(K) =  UserRate(Mopt(K),K);
    end
    Pi = zeros(1,KgOpt);
    for K=1:KgOpt
        [Pi(K),~] = MGmm_SD_Queue(K,KgOpt,(loading/100)*lambdaS,1,Ropt);
    end
    M_iterate(1)=round(sum(Mopt.*Pi));    
    if all(M_iterate == M_iterate(1))
        Converged = true;
        AvgEEperBS = sum(EEopt.*Pi);
        AvgURperBS = sum(Ropt.*Pi);
        disp(['Mopt_perUS = ' num2str(Mopt)])
        Outage = Pi(end)
    else
        M_iterate = ones(1,19).*min(M_iterate);
        disp(['M_it = ' num2str(M_iterate)]);
    end
    
end
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
end

