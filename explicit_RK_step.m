%This function computes the value of X at the next time step
%for any arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
%              have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
%           BT_struct.A: matrix of a_{ij} values
%           BT_struct.B: vector of b_i values
%           BT_struct.C: vector of c_i values
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
%    formula depends on the integration method used
%num_evals: A count of the number of times that you called
%           rate_func_in when computing the next step
function [XB, fun_count] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    % unpack the Butcher tableau
    A = BT_struct.A;
    B = BT_struct.B;
    C = BT_struct.C;
    
    % set up a vector to store K values
    K = zeros(length(C), length(XA));
    fun_count = 0;

    % loop through tableau values
    for i = 1:length(C)
        K(i,:) = rate_func_in(t+C(i)*h, XA + h*(A(i,:)*K));
        % add one to the function counter
        fun_count = fun_count + 1;
    end
    
    % calculate the next X value
    XB = XA + h*(B(1,:)*K);
end