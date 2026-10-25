% Function to plot local truncation error of the two embedded RK method
% outputs, the difference between these two outputs, 
%INPUTS:
% func: a function handle that defines a rate function to approximate
% sol: a function handle that defines the solution of func
% t_span: a two-value list containing the first and final values desired of
%         the relevant time span
% num_iter: the number of iterations/datapoints desired
%OUTPUTS:
% None
function plot_local_truncation(method, func, sol, t_ref, num_iter)
    % set optional inputs
    arguments
        method (1,:) string = "Heun Euler"
        func (1,:) function_handle = @rate_func01;
        sol (1,:) function_handle = @solution01;
        t_ref (1,:) double = 0.439
        num_iter (1,:) double = 100
    end

    % set methods
    data_color = [0 0 0; 0.5 0.5 1; 0.5 1 0.5; 1 0.5 0.5];
    data_titles = ["Ground truth $|X(t+h)-X(t)|$","$X_{diff}$ ($|XB1 - XB2|$)",...
        "$XB1$ Error ($|XB1-X_{true}|$)",...
        "$XB2$ Error ($|XB2-X_{true}|$)"];

    % set initial x guess and timestep size
    XA = sol(t_ref);
    h_list = logspace(-4,1,num_iter);

    % set up the Butcher tableau corresponding to the current method
    % loop, as well as empty lists for storing data this loop
    BT_struct = init_BT(method);
    X_diff = zeros(num_iter,1);
    true_diff = zeros(num_iter,1);
    XB1_error = zeros(num_iter,1);
    XB2_error = zeros(num_iter,1);

    % loop through different step sizes and calculate the truncation
    % errors for each one
    for i = 1:num_iter
        h = h_list(i);
        [XB1,XB2] = embedded_RK_step(func,t_ref,XA,h,BT_struct);

        X_true = sol(t_ref+h);
        true_diff(i) = norm(X_true - XA);
        X_diff(i) = norm(XB1 - XB2);
        XB1_error(i) = norm(XB1 - X_true);
        XB2_error(i) =  norm(XB2 - X_true);
    end

    % store all the data points in one list
    y_data = [true_diff X_diff XB1_error XB2_error];
    
    % set up unique filters to best capture the linear portions for each
    % graph
    x_true_filter = h_list <= 1e-3 & h_list >= 1e-4;
    x_diff_filter = h_list <= 1e-1 & h_list >= 1e-2;
    xb1_filter = h_list <= 1e-1 & h_list >= 1e-2;
    xb2_filter = h_list <= 1e-1 & h_list >= 1e-2;
    h_filter = [x_true_filter' x_diff_filter' xb1_filter' xb2_filter'];

    clf
    figure(1)

    for j = 1:4
        % loop through the list
        error = y_data(:,j);
        log_h_list = log10(h_list(h_filter(:,j)))';

        % p-value and slope calculations
        log_error = log10(error(h_filter(:,j)));
        p = polyfit(log_h_list, log_error, 1);
        disp(p(1))

        loglog(h_list, error, '.', 'markersize', 15, 'Color', data_color(j,:), ...
            'Displayname',sprintf("%s",data_titles(j)));
        hold on; % call hold on after calling loglog
        grid on;
    end

    % plot styling
    set(gca,"TickLabelInterpreter",'latex',fontsize=14);
    title(sprintf('Embedded RK (%s Method): Local Truncation Error vs Timestep',method), 'Interpreter', 'Latex', 'FontSize',18);
    legend('Location','southeast', 'Interpreter', 'Latex','FontSize',12);
    xlabel("Step size $h$ (-)", 'Interpreter', 'Latex','FontSize',16);
    ylabel("Local Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
    fontsize(scale=1.5)
    axis([h_list(1) h_list(end) 1e-17 1e6])    
end