% Optimization Driver 
clear; close all; clc;
rng = ('shuffle');

% following: 
tic
nvars = 9; 
Mutation_rate = 0.95;
Population_num = 300; 
Max_gen = 600;
% If not inside a slurm job, use 4 processors
% if isempty(getenv('SLURM_NTASKS'))
%     NP = 4;
% 
% else 
%     setupUmichClusters
%     NP = str2double(getenv('SLURM_NTASKS'));
%     poolobj = parpool('current', NP)
% 
% end

%  poolobj = parpool(15)

% genetic algorithm optimization options
gaoptions = optimoptions('ga','MaxGenerations',Max_gen,'Display','iter');
gaoptions = optimoptions(gaoptions,'UseParallel',true);
gaoptions = optimoptions(gaoptions,'PopulationSize',Population_num);
gaoptions = optimoptions(gaoptions,'FunctionTolerance',1e-6);
gaoptions = optimoptions(gaoptions,'MutationFcn', {@mutationadaptfeasible, Mutation_rate});
gaoptions = optimoptions(gaoptions,'MaxStallGenerations',1000);
gaoptions = optimoptions(gaoptions,'OutputFcn',@GA_DISP)
%% Simulation time span in Second:
% K_T = adjvar(1) * par(18); 
% kf = adjvar(2) * kf;
% k3 = adjvar(3) * k3;
% V0 = adjvar(4)* 0.4;
% n = adjvar(5); 
% f_myofibril = adjvar(6); % Percent Myofibril in Muscle from Palmer etal (Mol Cell Biochem. 2004 Aug;263(1-2):73-80)
% I have divided the Fpassive by a factor of 2, to hav almost the same
% range of passive force we have in the data. 

var = [2 1 1 1 2 1 1 1 1 1];

lb = var - var * 0.8;
ub = var + var * 20;
ConstraintFunction = [];
objectiveGA=@(x)objective_fun_multi(x);

% lb(7:8) = 1;
% ub(7:8) = 1;

adjvar = ga(objectiveGA,nvars,[],[],[],[],lb,ub,ConstraintFunction,gaoptions)

save run



% ATP2mM =[6.0000    1.0962
%     1.0000   27.4392
%     3.0000    4.3771
%     5.0000    1.8673
%     0.5000   43.5430
%     4.0000    2.4951
%     2.0000    7.6647
%          0   65.0903];
% 
% ATP8mM = [6	2.432467797
% 1	42.88019347
% 3	12.39339873
% 5	3.591327119
% 0.5	61.01338797
% 4	6.290181818
% 2	19.02498679
% 0	67.48992526];
% 
% force_ATP2mM_sorted = sort(ATP2mM(:,2),'descend')
% force_ATP8mM_sorted = sort(ATP8mM(:,2),'descend')
