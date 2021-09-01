%% Script to simulate Force-velocity data of Palmer etal J Mol Cell Cardiol. 2013 Apr;57:23-31
clear; close all;
global SLset timrel
SLset  = 2.2; % Set sarcomere length, Units: um
timrel = 0.1; % Time of sarcomere release, Units: ms
AxisFontSize = 18; LabelFontSize = 18;

% Set temperature and initial SL
TmpC = 21; SL0 = 2.2; % um

% Set metabolite concentrations, 
MgATP_All = [2 8]; MgADP = 0.0; Pi = 0.0; % Experimental conditions from Palmer etal J Mol Cell Cardiol. 2013 Apr;57:23-31
for k =1:length(MgATP_All)
    
MgATP = MgATP_All(k);
plot_color= {'b','k'};
init = [zeros(1,9),SL0]; % Initial conditions for the model
% running to get the F_Max
tspan = [0:0.0005:0.1]; nn = length(tspan);
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
[T, Y] = ode15s(@Model_XB,tspan,init,options,TmpC,MgATP,Pi,MgADP,0);
SL = Y(:,10); 
   
for i = 1:length(T)
        [~,~,Ftotal_0(i,k)] = Model_XB(T(i),Y(i,:),TmpC,MgATP,Pi,MgADP,0);
end

Fmax = Ftotal_0(end,k)
tspan = [0:0.0005:0.4]; nn = length(tspan);

% Fmax = 28.55; % Maximum Tension when sarcomere stretched to 2.2 um at 17 C
F_load = [0.0:0.05:1.0]*Fmax; % Afterloads agains which the sarcomere contracts
mm = length(F_load); vm = zeros(1,mm);
dSL = zeros(nn,mm); Ftotal = zeros(nn,mm);

tic;
for j = 1:mm
    [T, Y] = ode15s(@Model_XB,tspan,init,options,TmpC,MgATP,Pi,MgADP,F_load(j));
    SL = Y(:,10); 
   
    for i = 1:length(T)
        [~, dSL(i,j),Ftotal(i,j)] = Model_XB(T(i),Y(i,:),TmpC,MgATP,Pi,MgADP,F_load(j));
    end
        ts = 321; X = ts-10:ts+10; % time-points at which the average slope of sliding velocity is calculated i.e. between 40 to 50 ms after release
        vm(1,j) = mean(dSL(X,j)./SL0);
        figure(2); set(figure(2),'Units','inches','Position',[11 0.5 5 7]);
        subplot(211), plot(T,SL,plot_color{k},'linewidth',1.5); ylim([1.4 2.35]); hold on; ylabel('SL (\mu m)','fontsize',LabelFontSize)
        title(sprintf('Tmp = %d C',TmpC),'fontsize',14); xlabel('Time (sec)','fontsize',LabelFontSize), set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
        subplot(212), plot(T,Ftotal(:,j),plot_color{k},'linewidth',1.5); hold on; 
        ylabel('Force (mN mm^{-2})','fontsize',LabelFontSize)
        xlabel('Time (sec)','fontsize',LabelFontSize)
        set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
clear T Y SL
end
toc;
%% Plots Figure 4 of the article
figure(1); 
% clf; 
set(figure(1),'Units','inches','Position',[0.5 0.5 10 5]);
subplot('position',[0.1 0.15 0.35 0.8])
% plot(F_load/Fmax,-(vm),'linewidth',2.5,'color',plot_color{k});
plot(F_load,-(vm),'linewidth',2.5,'color',plot_color{k}); 
% xlim([0 1]); 
hold on;

ylabel('Velocity (SL s^{-1})','fontsize',LabelFontSize)
xlabel('T/Tmax','fontsize',LabelFontSize)
set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
% text(-0.15, 1.40,'A.','fontsize',22)

% subplot('position',[0.55 0.15 0.35 0.8]), plot(F_load/Fmax,-(vm).*F_load/Fmax,'linewidth',2.5,'color',plot_color{k}); 
subplot('position',[0.55 0.15 0.35 0.8]), plot(F_load,-(vm).*F_load,'linewidth',2.5,'color',plot_color{k}); 

hold on; 
% ylim([0.0 0.3]); 
hold on; 
xlabel('T/Tmax','fontsize',LabelFontSize), ylabel('Normalised Power','fontsize',LabelFontSize)
% text(-0.15, 0.29,'B.','fontsize',22)
set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);

%% Plots Data from Palmer etal J Mol Cell Cardiol. 2013 Apr;57:23-31
data_x = [0.0391    0.0797    0.1300    0.1580    0.1970    0.2480    0.3010    0.3490    0.3960    0.4990    0.5990    0.7000    0.8000    0.9000];
data_y = [1.3000    1.2000    1.0000    0.9600    0.9000    0.7900    0.7000    0.5900    0.5200    0.4000    0.2900    0.1900    0.0910    0.0180];
figure(1);
% subplot('position',[0.1 0.15 0.35 0.8]),plot(data_x,data_y,'o','color',plot_color{k},'MarkerEdgeColor',plot_color{k}, ...
%     'MarkerFaceColor','w', 'Markersize',10,'linewidth',2); xlim([0 1]); hold on;
ylabel('Velocity (L s^{-1})','fontsize',LabelFontSize)
xlabel('T/Tmax','fontsize',LabelFontSize)
set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize); 
% l = legend('Model','Data'); 
l = legend('ATP = 2 mM','ATP = 8 mM'); 

set(l,'interpreter','latex','Fontsize',LabelFontSize);
% ylim([0 1.45]); 
hold on;
% subplot('position',[0.55 0.15 0.35 0.8]), plot(data_x,data_x.*data_y,'o','color',plot_color{k},'MarkerEdgeColor',plot_color{k}, ...
%     'MarkerFaceColor','w', 'Markersize',10,'linewidth',2); hold on; 
% ylim([0.0 0.3]); 
hold on; 
xlabel('T/Tmax','fontsize',LabelFontSize)
ylabel('Normalised Power','fontsize',LabelFontSize);  
% l = legend('Model','Data');
l = legend('ATP = 2 mM','ATP = 8 mM'); 

set(l,'interpreter','latex','Fontsize',LabelFontSize);
set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
end