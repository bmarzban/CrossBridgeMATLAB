%% Script to simulate Force-velocity data of Palmer etal J Mol Cell Cardiol. 2013 Apr;57:23-31
%  clear; close all;
global SLset timrel
SLset  = 2.2; % Set sarcomere length, Units: um
timrel = 0.1; % Time of sarcomere release, Units: ms
AxisFontSize = 18; LabelFontSize = 18;

% Set temperature and initial SL
TmpC = 21; SL0 = 2.2; % um

% Set metabolite concentrations, 
MgATP = 2.0; MgADP = 0.0; Pi = 0.0; % Experimental conditions from Palmer etal J Mol Cell Cardiol. 2013 Apr;57:23-31

init = [zeros(1,9),SL0]; % Initial conditions for the model
Velocity = [0.0:-0.5:-4.5]; % Afterloads agains which the sarcomere contracts
% Velocity = -4.4;
mm = length(Velocity); vm = zeros(1,mm);
tspan = [0:0.0005:0.7]; nn = length(tspan);
dSL = zeros(nn,mm); Ftotal = zeros(nn,mm);

tic;
for j = 1:mm
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    [T, Y] = ode15s(@Model_XB,tspan,init,options,TmpC,MgATP,Pi,MgADP,Velocity(j));
    SL = Y(:,10); 
       
    for i = 1:length(T)
        [~, dSL(i,j),Ftotal(i,j)] = Model_XB(T(i),Y(i,:),TmpC,MgATP,Pi,MgADP,Velocity(j));
    end
    figure(1)
    hold on
    plot(T, Ftotal(:,j))
    xlabel('time (s)')
    ylabel('Force (kPa)')
    figure(2)
    hold on
    plot(T,SL)
    xlabel('time (s)')
    ylabel('SL (um)')
%    legend(num2str(Velocity(j)))
% clear T Y SL
%  Force (j) = mean(Ftotal(find(SL < 2.01 & SL>1.99),j))
%  if j==1;
%  Force (1) = max(Ftotal(:,j));
Force(j) = Ftotal(end,j);
%  end
end
figure(1)
legend( '0'       , '-1 um/s'  ,        '-2 um/s'    ,    '-3 um/s'    ,      '-4um/s',    '-5 um/s'    ,      '-6um/s')
figure(2)
legend( '0'       , '-1 um/s'  ,        '-2 um/s'    ,    '-3 um/s'    ,      '-4um/s',    '-5 um/s'    ,      '-6um/s')
% figure(3)
% hold on
% plot(Force,abs(Velocity).*Force)
% xlabel('Force (kPa)')
% ylabel('Power (Watts)')
% legend('ATP = 2mM','ATP = 8mM')

figure(4)
hold on
plot(Force,abs(Velocity))
xlabel('Force (kPa)')
ylabel('Velocity (um/s)')
legend('ATP = 2mM','ATP = 8mM')

% ylim([0 5])
xlim([0.0 100])