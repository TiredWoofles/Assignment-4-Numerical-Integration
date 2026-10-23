%Runs numerical integration arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
%              have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
%           BT_struct.A: matrix of a_{ij} values
%           BT_struct.B: vector of b_i values
%           BT_struct.C: vector of c_i values
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, fun_count] = explicit_RK_fixed_step_integration ...
        (rate_func_in,tspan,X0,h_ref,BT_struct)
    % extract t0 and tf
    t0 = tspan(1);
    tf = tspan(2);
    
    % find smallest N that's smaller than h_ref
    N = 1;
    while (tf-t0)/N >= h_ref
        N = N+1;
    end 

    % calculate the largest h value that will be an even fraction of the
    % time span
    h_avg = (tf-t0)/N;
    
    % initiate parameters for running the step function
    t_list = linspace(t0,tf,N+1)';
    X_list = zeros(N+1,length(X0));
    % transpose X0 for multidimensional inputs
    X_list(1,:) = X0(:)';
    fun_count = 0;
    
    % loop through number of values N
    for j = 1:N
        t = t_list(j);
        XA = X_list(j,:);
        
        % solve the step function
        [XB, iter] = explicit_RK_step(rate_func_in,t,XA,h_avg,BT_struct);
    
        % set the next list value as the solution, and add one to the
        % function counter
        X_list(j+1,:) = XB;
        fun_count = fun_count + iter;
    end
end