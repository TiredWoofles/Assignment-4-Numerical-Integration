%This function computes the value of X at the next time step
%using the implicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB,num_evals] = implicit_midpoint_step(rate_func_in,t,XA,h)
    num_evals = 0;
    G = @(XB) count_evals(XB);
    function G_value = count_evals(XB)
        num_evals = num_evals +1;
        G_value = XA +h*rate_func_in(t + h/2, (XA +XB)/2)-XB;
    end
    solver_params.numerical_diff = 1;
    [XB,exit_flag] = multi_newton_solver(G,XA,solver_params);
end
