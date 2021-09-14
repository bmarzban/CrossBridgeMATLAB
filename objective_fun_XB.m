function err = objective_fun_XB(adjvar,MgATP)
% Set temperature and initial SL
TmpC = 21; SL0 = 2.2; % um

% Set metabolite concentrations, 
% MgATP = 8.0; 
MgADP = 0.0;
Pi = 0.0; % Experimental conditions from Palmer etal J Mol Cell Cardiol. 2013 Apr;57:23-31

init = [zeros(1,9),SL0]; % Initial conditions for the model
Velocity = [0, -0.5,-1, -2 ,-3,-4,-5]; % Afterloads agains which the sarcomere contracts
% % Velocity = [0.0:-0.5:-6]; % Afterloads agains which the sarcomere contracts


force_ATP2mM_sorted = [65.0903   43.5430   27.4392    7.6647    4.3771    2.4951    1.8673   ];
force_ATP8mM_sorted = [67.4899   61.0134   42.8802   19.0250   12.3934    6.2902    3.5913   ];
if MgATP == 2
    force_target = force_ATP2mM_sorted;
elseif MgATP == 8
    force_target = force_ATP8mM_sorted;
end
mm = length(Velocity); vm = zeros(1,mm);
tspan = [0:0.0005:0.25]; nn = length(tspan);
dSL = zeros(nn,mm); Ftotal = zeros(nn,mm);
for j = 1:mm
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    [T, Y] = ode15s(@Model_XB,tspan,init,options,TmpC,MgATP,Pi,MgADP,Velocity(j),adjvar);
%     SL = Y(:,10); 
       
    for i = 1:length(T)
        [~, dSL(i,j),Ftotal(i,j)] = Model_XB(T(i),Y(i,:),TmpC,MgATP,Pi,MgADP,Velocity(j),adjvar);
    end
   
Force(j) = Ftotal(end,j);
%  end
end
 err = sumsqr(Force - force_target);
