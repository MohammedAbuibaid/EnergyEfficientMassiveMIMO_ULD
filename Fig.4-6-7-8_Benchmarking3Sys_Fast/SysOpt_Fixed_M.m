function [AvgEEperBS, AvgURperBS] = SysOpt_Fixed_M(PLO_Network,PLI_Network,KgOpt,MgOpt,Pc,PmaxPA,loading,lambdaS)

UR_per_BS= zeros(19,1);
EE_per_BS=  zeros(19,1);

activityRef=ones(1,19);
Converged=false; 
RoundNum = 1;
% disp(['activityRef = ' num2str(activityRef)])

while (~Converged)
    %disp(['$$$$$$$$$$$$$$$$$$$$$   Round Num   ' num2str(RoundNum) '   $$$$$$$$$$$$$$$$$$$$$']);
    activityIterate=activityRef;
    %disp(['activityIte = ' num2str(activityIterate)])
    
    for BSj= 1:19
         %disp('---------------------------------------------------------------------------------------------------')
         %disp(['For BSj = ' num2str(BSj)])
               
        % Calculating SIR according to Activity Factor
        ActivityTemp=activityIterate;
        ActivityTemp(BSj)=[];
        
        EE = zeros(1,KgOpt);
        UserRate = zeros(1,KgOpt);
        Mc=MgOpt; Md = ActivityTemp.*MgOpt.*ones(1,18);
        p = Pc/MgOpt;
        for K=1:KgOpt
            [EE(K),UserRate(K)] = EE_R_Ptot_PA(PLO_Network{BSj},PLI_Network{BSj},KgOpt,K,Mc,Md,p,PmaxPA,[]);
        end
        
        Pi = zeros(1,KgOpt);
        for K=1:KgOpt
            [Pi(K),P0] = MGmm_SD_Queue(K,KgOpt,(loading/100)*lambdaS,1,UserRate);
        end
        % Probability of having non-zero users in the BS(ww)
        
        activityIterate(BSj)= (1-P0);
        %disp(['activityIterate of BSj  = ' num2str(1-P0)])    
        %disp(['activityIte = ' num2str(activityIterate)])

        EE_per_BS(BSj) = sum(EE.*Pi);
        UR_per_BS(BSj)= sum(UserRate.*Pi);
    end 
    
    epsilon= max(abs(activityRef-activityIterate));
    % disp(['epsilon = ' num2str(epsilon)])

    if epsilon <10^-60
        Converged=true; % All BSs are Converged
%          disp(['$$$$$$$$$$$$$$$$$$$$$   System Converged After ' num2str(RoundNum) ' $$$$$$$$$$$$$$$$$$$$$'])
%          disp(['epsilon = ' num2str(epsilon)])
%          disp(['activityIte = ' num2str(activityIterate)])
               
    else
        
        RoundNum = RoundNum+1;
%          disp(['activityRef = ' num2str(activityRef) ])
        activityRef= activityIterate;
%          disp('Update activityRef to activityIterate')
%          disp(['activityRef = ' num2str(activityRef) ])
        
    end
    
end % While

AvgURperBS=mean(UR_per_BS);
AvgEEperBS = mean(EE_per_BS);
Outage = Pi(end);


end