%Runs numerical integration using forward Euler approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration

function [t_list,X_list,h_avg, num_evals] = implicit_midpoint_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
  %Initialise the lists
  t_list = [];
  X_list = [];
  num_evals = 0;
  
  %Find h_avg
  t0 = tspan(1);
  tf = tspan(2);

  N = (tf-t0)/h_ref; %current number of steps
  N = ceil(N); %Round to nearest whole number

  h_avg = (tf-t0)/N; %Actual step size
   
  X_current = X0;

  t_span = linspace(t0, tf, N); %Generate list of times
  for i = 1:length(t_span)
      t = t_span(i);
      t_list(i) = t; %Add to time vector
      [Xn, num_evals_n] = implicit_midpoint_step(rate_func_in, t, X_current, h_avg); %Evaluate new X
      num_evals = num_evals + num_evals_n;
      X_list(:, i) = Xn; %Add to X value vector
      X_current = Xn;
  end
end