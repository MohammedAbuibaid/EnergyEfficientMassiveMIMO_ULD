function [Mopt,AvgEEperBS, AvgURperBS, AvgPtotperBS, Mavg]=SysOpt_Adaptive_M_SingleCell(PLO_Network,PLI_Network,KgOpt,Mmin,MgOpt,Pc,PmaxPA,loading,lambdaS)

disp(['++++ Loading = ' num2str(loading) '% ++++'])
% disp('For BSj = 1');
Mref=MgOpt;
% disp(['Mref = ' num2str(Mref) ])
Converged=false;
M_iterate=Mref;

while(~Converged)
    
    % Finding the optimum antennas for each user state
    Md_temp= M_iterate*ones(1,18);
    EE = zeros(MgOpt,KgOpt);
    UserRate = zeros(MgOpt,KgOpt);
    Ptot = zeros(MgOpt,KgOpt);
    for K=1:KgOpt
        for M=Mmin:MgOpt
            Mc = M; Md = Md_temp;
            p = Pc/M;
            [EE(M,K),UserRate(M,K),Ptot(M,K)] = EE_R_Ptot_PA(PLO_Network,PLI_Network,KgOpt,K,Mc,Md,p,PmaxPA,[]);
        end
    end
    [EEopt, Mopt]=max(EE);
    
    
    Ropt = zeros(1,KgOpt);
    Ptot_opt = zeros(1,KgOpt);
    for K=1:KgOpt
        Ropt(K) =  UserRate(Mopt(K),K);
        Ptot_opt(K) =  Ptot(Mopt(K),K);
    end
    
    
    Pi = zeros(1,KgOpt);
    for K=1:KgOpt
        [Pi(K),~] = MGmm_SD_Queue(K,KgOpt,(loading/100)*lambdaS,1,Ropt);
    end
    
    
    Mavg = ceil(sum(Mopt.*Pi));
    
        
    if  Mavg == M_iterate
        Converged = true;
        disp(['Mavg = ' num2str(Mavg) ])
        %disp(['Mopt(end) = ' num2str(Mopt(end))])
        AvgEEperBS = sum(EEopt.*Pi);
        AvgURperBS = sum(Ropt.*Pi);
        AvgPtotperBS = sum(Ptot_opt.*Pi);
        %disp(['Mopt = ' num2str(Mopt)])
        % Outage = Pi(end);
    else
        M_iterate = Mavg;
        %disp(['Mopt = ' num2str(Mopt)])
        %disp(['M_it = ' num2str(M_iterate)]);
        
    end
    
end
end

