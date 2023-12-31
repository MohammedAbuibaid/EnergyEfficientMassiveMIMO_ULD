clear all; close all; clc;
load('Res_C2B_Design.mat');
B_Ratio = Cen2Bou_Ratio(2,:);

%% KgOpt Vs. MgOpt
figure
hold on; 
set(gca,'fontsize',10)
[hAx,hLine1,hLine2] = plotyy(B_Ratio,KgOpt,B_Ratio,MgOpt);
% title('KgOpt Vs. MgOpt')
x = {'100/0' '90/10' '80/20' '70/30' '60/40' '50/50' '40/60' '30/70' '20/80' '10/90' '0/100'};
set(gca,'xticklabel',x.')
xlabel('Center/Boundary Users Ratio (%)')
ylabel(hAx(1),'K_{max}')
ylabel(hAx(2),'M_{max}')
hLine1.LineStyle = '--';
hLine2.LineStyle = '--';
hLine1.Marker = 'd';
hLine2.Marker = 'o';
hLine1.LineWidth = 1.2;
hLine2.LineWidth = 1.2;
set(gcf, 'Position', [500 500 900 500])


%% EE vs. User Rates Rc Rb Rt

Rb_Gain = 100*(Rb(11)-Rb(2))/Rb(2);
Rc_Gain = 100*(Rc(11)-Rc(2))/Rc(2);

disp(['Rb Gain = ' num2str(Rb_Gain)]);
disp(['Rc Gain = ' num2str(Rc_Gain)]);

%%
figure
hold on; 
set(gca,'fontsize',10)
% [hAx,hLine1,hLine2] = plotyy(B_Ratio,EE_MK_design/1e6,[B_Ratio',B_Ratio',B_Ratio'],[Rb',Rc',Rt']./1e6);
[hAx,hLine1,hLine2] = plotyy(B_Ratio,EE_MK_design/1e6,[B_Ratio',B_Ratio'],[Rb',Rc']./1e6);
% title('BS Energy Efficiency Vs. User Rates ( Rc Rb ) ')
x = {'100/0' '90/10' '80/20' '70/30' '60/40' '50/50' '40/60' '30/70' '20/80' '10/90' '0/100'};
set(gca,'xticklabel',x.')
xlabel('Center/Boundary Users Ratio (%)')
ylabel(hAx(1),'Energy Efficiency (Mbit/Joule)')
ylabel(hAx(2),'User Rate (Mbps)')
hLine1.LineStyle = '--';
hLine1.LineWidth = 1.2;
hLine1.Color = 'b';
hLine1.Marker = 's';

hLine2(1).LineStyle = '-.';
hLine2(2).LineStyle = '-.';
%hLine2(3).LineStyle = '-.';
hLine2(1).LineWidth = 1.2;
hLine2(2).LineWidth = 1.2;
%hLine2(3).LineWidth = 1.2;
hLine2(1).Color = 'r';
hLine2(2).Color = 'g';
%hLine2(3).Color = 'k';
hLine2(1).Marker = 'x';
hLine2(2).Marker = 'o';
%legend('EE','Rb','Rc','Rt','location','northwest')
legend('EE','Rb','Rc','location','best')
set(gcf, 'Position', [500 500 900 500])



% %% Weighted Energy Efficiency Vs. Weighted User Rate
% figure
% hold on; grid on;
% set(gca,'fontsize',10)
% [hAx,hLine1,hLine2] = plotyy(B_Ratio,WA_EE/1e6,B_Ratio,WA_R/1e6);
% title('BS Energy Efficiency Vs. User Rate')
% x = {'100/0' '90/10' '80/20' '70/30' '60/40' '50/50' '40/60' '30/70' '20/80' '10/90' '0/100'};
% set(gca,'xticklabel',x.')
% xlabel('Center/Boundary Users Ratio %')
% ylabel(hAx(1),'Weighted EE [Mbit/Joule]')
% ylabel(hAx(2),'Weighted User Rate [Mbps]')
% hLine1.LineStyle = '--';
% hLine2.LineStyle = '--';
% hLine1.Marker = 'd';
% hLine2.Marker = 'o';
% hLine1.LineWidth = 1.2;
% hLine2.LineWidth = 1.2;
%%