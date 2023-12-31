function [Mopt,AvgEEperBS, AvgURperBS]=SysOpt_Adaptive_M(PLO_Network,PLI_Network,KgOpt,Mmin,MgOpt,Pc,PmaxPA,loading,lambdaS)

EE_per_BS =zeros(19,1);
UR_per_BS =zeros(19,1);

Mref=MgOpt*ones(1,19);
disp(['Mref = ' num2str(Mref) ])
Converged=false;
M_iterate=Mref;
RoundNum = 1;

while (~Converged)
    
     disp(['$$$$$$$$$$$$$$$$$$$$$   Round Num   ' num2str(RoundNum) '   $$$$$$$$$$$$$$$$$$$$$']);
     disp(['M_it = ' num2str(M_iterate) ])
     
    for BSj=1:19
        disp('---------------------------------------------------------------------------------------------------')
        disp(['For BSj = ' num2str(BSj)])
        Md_temp= M_iterate;
        Md_temp(BSj)= [];
        
        % Finding the optimum antennas for each user state
        EE = zeros(MgOpt,KgOpt);
        UserRate = zeros(MgOpt,KgOpt);
        for K=1:KgOpt
            for M=Mmin:MgOpt;
                Mc = M; Md = Md_temp;
                p = Pc/M;
                [EE(M,K),UserRate(M,K)] = EE_R_Ptot_PA(PLO_Network{BSj},PLI_Network{BSj},KgOpt,K,Mc,Md,p,PmaxPA,[]);
            end
        end
        [EEopt, Mopt]=max(EE);
        disp(['Mopt_perUS = ' num2str(Mopt)])
        Ropt = zeros(1,KgOpt);
        for K=1:KgOpt
            Ropt(K) =  UserRate(Mopt(K),K);
        end
        Pi = zeros(1,KgOpt);
        for K=1:KgOpt
            [Pi(K),~] = MGmm_SD_Queue(K,KgOpt,(loading/100)*lambdaS,1,Ropt);
        end
        M_iterate(BSj)=round(sum(Mopt.*Pi));
        disp(['Weighted-Avg Mopt =  ' num2str(round(sum(Mopt.*Pi)))]);
        disp(['M_it = ' num2str(M_iterate)]);
        EE_per_BS(BSj) = sum(EEopt.*Pi);
        UR_per_BS(BSj)= sum(Ropt.*Pi);
        
    end
    
    if (Mref-M_iterate)==0
        Converged=true;
         disp(['$$$$$$$$$$$$$$$$$$$$$   System Converged After ' num2str(RoundNum) ' $$$$$$$$$$$$$$$$$$$$$'])
         disp(['Mopt_perBS = ' num2str(M_iterate) ])
    else
        RoundNum = RoundNum+1;
        disp(['M_ref = ' num2str(Mref) ])
        Mref=M_iterate;
        disp('Update Mref to M_iterate')
        disp(['M_ref = ' num2str(Mref) ])
    end
    
end

AvgEEperBS = mean(EE_per_BS);
AvgURperBS=mean(UR_per_BS);
Outage = Pi(end)
end

