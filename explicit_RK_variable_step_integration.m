%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, fun_count] = explicit_RK_variable_step_integration ...
(func,tspan,X0,h_ref,BT_struct,p,error_desired)
    % initialize parameters
    t0 = tspan(1);
    tf = tspan(2);
    t = t0;
    XA = X0;
    h_avg = h_ref;
    fun_count = 0;

    % create empty lists for saving timesteps and X values
    t_list = t;
    X_list = XA;
    
    % loop through number of values N
    while t < tf
        XA = X_list(end,:);

        % solve the step function
        [XB, iter, h_next, redo] = explicit_RK_variable_step(func,t,XA,h_avg, ...
            BT_struct,p,error_desired,growth_limit);
    
        if ~redo
            % set the next list value as the solution, and add one to the
            % function counter
            X_list(end+1,:) = XB';
            if t + h_next > tf
                h_next = tf - t;
            end
            t = t+h_next;
            t_list(end+1) = t;
        end

        fun_count = fun_count + iter;
        h_next = h_avg;
    end
end