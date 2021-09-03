function [state,options,optchanged] = GA_DISP(options,state,flag)

optchanged = false;
 
switch flag
    
    case 'iter'
        ibest = state.Best(end);
        ibest = find(state.Score == ibest,1,'last');
        bestx = state.Population(ibest,:);
        
%         load inputs inputs
%         PlotObjectiveMechanoEnergetics(bestx,inputs)
%         filename_mat1 = strcat(pwd,'/bestGroupsFits',num2str(inputs(1)),'.mat');
% 
%         load (filename_mat1,'bestGroupsFits')
%         bestGroupsFits = [bestGroupsFits; bestx];                                % Read Previous Results, Append New Value
%         save(filename_mat1, 'bestGroupsFits')       
       
        disp(['Best GA: ', num2str(bestx)]);
       
    otherwise
        disp(['intial or done =) ']);
end
