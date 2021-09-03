function err_total = objective_fun_multi(adjvar)
try
MgATP = 2;
err2 = objective_fun_XB(adjvar,MgATP);
MgATP = 8;
err8 = objective_fun_XB(adjvar,MgATP);
catch 
    err2 = 10e5;
    err8 = 5e5;
end
% ramp_err = objective_function_step_input(adjvar,F_exp_spline);
% err_total = err2 + err8 + ramp_err;
err_total = err2 + err8;