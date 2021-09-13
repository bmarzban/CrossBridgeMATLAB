%% Script to simulate Force-velocity data of Palmer etal J Mol Cell Cardiol. 2013 Apr;57:23-31
clear; close all;
AxisFontSize = 18; LabelFontSize = 18;
% Set temperature and initial SL
TmpC = 21; SL0 = 2.2; % um

% Set metabolite concentrations, 
MgATP_all = [2 8]; 
MgADP = 0.0; Pi = 0.0; % Experimental conditions from Palmer etal J Mol Cell Cardiol. 2013 Apr;57:23-31
for k = 1:length(MgATP_all)
    MgATP= MgATP_all(k);
init = [zeros(1,9),SL0]; % Initial conditions for the model
% Velocity = [0.0:-0.5:-6]; % Afterloads agains which the sarcomere contracts
Velocity = [0, -0.5,-1, -2 ,-3,-4,-5]
% Velocity = -4.4;
mm = length(Velocity); vm = zeros(1,mm);
tspan = [0:0.0005:0.7]; nn = length(tspan);

dSL = zeros(nn,mm); Ftotal = zeros(nn,mm);
adjvar = [2 1 1 1 2 1]
adjvar = [1.2288      0.9233     0.92438     0.96558      3.1201     0.78526]
adjvar = [1.6         0.8     0.94184         0.8         2.6     0.81503]
adjvar = [ 1.0456       1.096     0.52951     0.59008      2.1276     0.60655]
adjvar =[1.0001     0.50012     0.63637     0.70026      1.9596     0.70001]
adjvar =[1.2027      1.6527      1.0415     0.67868      2.4079     0.88576]
adjvar=[1.6         0.8     0.87832         0.8      2.0877         0.8]
adjvar=[1.6         0.8     0.87832         0.8      2.0877         0.8]
% adjvar=[0.50183     0.43513     0.61968     0.72101      1.1179     0.70746]
% adjvar = [0.51537     0.79514     0.48867     0.60496      1.1575     0.60655]
% adjvar = [0.5         0.4      0.4204         0.6      1.4422         0.6]

%k3stiff k2 and k1
adjvar=[1.21923      3.99625     0.226564      1.35106     0.839729     0.609746      13.8377      10.2223 1 1 1]
adjvar = [10.9139     0.801401      1.44883     0.815298      1.92399     0.823438          0.8      1.15156      10.9985      10.9985          0.8]
% adjvar = [1.0349      3.9881     0.23607      1.7181     0.73584     0.60303      5.9946      5.9907      5.9826     0.61526     0.95287]
% adjvar =[1.0506      3.9959     0.23607      1.6713     0.73584     0.60694      5.9985      5.9907      5.9983     0.60744      2.9529]
% adjvar = [1.61066     0.817956      5.41663      5.82525      21.9968     0.805524      10.9774      10.494 1 1 1]
% adjvar = [ 1.6          0.8      5.36597      5.82428           22          0.8           11           11 1 1 1]
 adjvar = [13.1632          0.8       1.6031          0.8      1.91318     0.902125          0.8      0.89256           11           11          0.8]
adjvar =[1.6          0.8      5.36591      5.82473           22          0.8           11           11 1 1 1]
adjvar =[1.60274      1.67384      6.84905      2.36059      2.54046     0.946731      10.9995      10.9986     0.801933]
tic;
for j = 1:mm
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',1e-3);
    [T, Y] = ode15s(@Model_XB,tspan,init,options,TmpC,MgATP,Pi,MgADP,Velocity(j),adjvar);
    SL = Y(:,10); 
       
    for i = 1:length(T)
        [~, dSL(i,j),Ftotal(i,j)] = Model_XB(T(i),Y(i,:),TmpC,MgATP,Pi,MgADP,Velocity(j),adjvar);
    end
    figure(1)
    hold on
    plot(T, Ftotal(:,j))
    xlabel('time (s)')
    ylabel('Force (kPa)')
%     figure(2)
%     hold on
%     plot(T,SL)
%     xlabel('time (s)')
%     ylabel('SL (um)')
%    legend(num2str(Velocity(j)))
% clear T Y SL
%  Force (j) = mean(Ftotal(find(SL < 2.01 & SL>1.99),j))
%  if j==1;
%  Force (1) = max(Ftotal(:,j));
Force(j) = Ftotal(end,j);

%  end

end
% 
% force_ATP2mM_sorted = [65.0903   43.5430   27.4392    7.6647    4.3771    2.4951    1.8673    1.0962];
% force_ATP8mM_sorted = [67.4899   61.0134   42.8802   19.0250   12.3934    6.2902    3.5913    2.4325];
force_ATP2mM_sorted = [65.0903   43.5430   27.4392    7.6647    4.3771    2.4951    1.8673   ];
force_ATP8mM_sorted = [67.4899   61.0134   42.8802   19.0250   12.3934    6.2902    3.5913   ];

if MgATP == 2
    force_target = force_ATP2mM_sorted;
elseif MgATP == 8
    force_target = force_ATP8mM_sorted;
end
 err = sumsqr(Force - force_target)

figure(1)
legend( '0'       , '-1 um/s'  ,        '-2 um/s'    ,    '-3 um/s'    ,      '-4um/s',    '-5 um/s'    ,      '-6um/s')
% figure(2)
% legend( '0'       , '-1 um/s'  ,        '-2 um/s'    ,    '-3 um/s'    ,      '-4um/s',    '-5 um/s'    ,      '-6um/s')
% % figure(3)
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
scatter(force_target,abs(Velocity))

% ylim([0 5])
xlim([0.0 100])
end
legend('ATP = 2mM','ATP = 2mM','ATP = 8mM','ATP = 8mM')
