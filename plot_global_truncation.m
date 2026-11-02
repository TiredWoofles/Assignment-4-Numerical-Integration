% Function to plot global truncation error against step size or function
% count for various Runge-Kutta methods.
%INPUTS (all optional):
% var_test: a string that calls function count or step size calculations;
%           defaults to step size
% func: a function handle that defines a rate function to approximate
% sol: a function handle that defines the solution of func
% t_span: a two-value list containing the first and final values desired of
%         the relevant time span
% num_iter: the number of iterations/datapoints desired
%OUTPUTS:
% None
function [error_list,fun_count_list] = plot_global_truncation(run, var_test, error_list, fun_count, func, sol, t_span, num_iter)
    arguments
        run (1,:) double = 1
        var_test (1,:) string = "Step Size"
        error_list (:,3) double = [0,0,0]
        fun_count (:,3) double = [0,0,0]
        func (1,:) function_handle = @gravity_wrapper
        sol (1,:) function_handle = @compute_planet_wrapper
        t_span (1,:) double = [0, 30]
        num_iter (1,:) double = 50
    end

    % set methods
    method = ["Van der Houwen's Third-Order","Kutta's Third-Order","Ralston's"];
    data_color = [0.65 .25 .10; 0.25 0 0.65; 0.25 0.65 0];
    fit_color = [1 .75 0.55; 0.65 0.35 1; 0.75 1 0];

    % set initial x guess, timestep size, a filter based on the timestep
    % size (for calculating p), and loop iterations
    X0 = sol(t_span(1));
    h_list = logspace(-5,1,num_iter);
    iter = 0;

    % set up the figure outside of the for-loop & if-statements
    figure(1)

    if run
        error_list = zeros(num_iter,3);
        fun_count_list = zeros(num_iter,3);
    end 

    % Loop through the three RK methods
    for j = 1:length(method)
        % set up the Butcher tableau corresponding to the curernt method: 
        % BT and filter boundaries
        [BT_struct,~,bounds,axis_range] = init_BT(method(j),var_test);
        if run
            % set up empty lists for storing data this loop
            error = zeros(num_iter,1);
            fun_count = zeros(num_iter,1);
    
            % loop through different step sizes and calculate the truncation
            % errors for each one
            for i = 1:num_iter
                h = h_list(i);
                [~, x_list,~, evals] = explicit_RK_fixed_step_integration ...
                    (func,t_span,X0,h,BT_struct);  
                x_true = sol(t_span(end));
                error(i) = norm(x_list(end,:).' - x_true);
                fun_count(i) = evals;
                
                % print iterations for troubleshooting
                iter = iter + 1
            end
            error_list(:,j) = error
        else
            error = error_list(:,j);
        end
    
        % Calculate p-value, plot a fit line, and plot data for a
        % comparison of error and function count.
        if var_test == "Fun Count"
            % set up a filter for the linear portion of the method
            filter = fun_count >= bounds(1) & fun_count <= bounds(2);
            log_error = log10(error(filter));
            log_evals = log10(fun_count(filter));

            % calculate the p-value (aka the slope) and fit line
            slope = polyfit(log_evals,log_error, 1);
            fit_line = 10.^polyval(slope, log10(fun_count(filter)));
            p = slope(1)

            % plot the data
            loglog(fun_count, fit_line, '-', 'Linewidth', 7, 'Color', fit_color(j,:), 'Displayname',sprintf("Fixed Step %s Fit",method(j)));
            hold on;
            grid on;
            loglog(fun_count, error, '.', 'MarkerSize',25, 'Color', data_color(j,:),'Displayname',sprintf("Fixed Step %s Data",method(j)));
            
            % plot styling
            set(gca,"TickLabelInterpreter",'latex');
            title('Global Truncation Error $\epsilon$ for RK  vs Efficiency', 'Interpreter', 'Latex', 'FontSize',18);
            legend('Location','northwest', 'Interpreter', 'Latex','FontSize',12);
            xlabel("Number of Rate Function Calls (-)", 'Interpreter', 'Latex','FontSize',16);
            ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
            legend('Location','northeast', 'Interpreter', 'Latex','FontSize',10);
            fontsize(scale=1.5)
            axis([10 6e6 18e-10 75])
    
        % Alternatively, calculate p-value, plot a fit line, and plot data 
        % for a comparison of error and step size.
        elseif var_test == "Step Size"
            % set up a filter for the linear portion of the method
            filter = h_list >= bounds(1) & h_list <= bounds(2);
            log_error = log10(error(filter));
            log_h_list = log10(h_list(filter));

            % calculate the p-value (aka the slope) and fit line
            slope = polyfit(log_h_list, log_error, 1);
            h_fit = h_list(filter);
            fit_line = 10.^polyval(slope, log10(h_fit));
            p = slope(1)
    
            % plot the data
            loglog(h_fit, fit_line, '-', 'Linewidth', 7, 'Color', fit_color(j,:), 'Displayname',sprintf("%s Fit",method(j)));
            hold on; % call hold on after calling loglog
            loglog(h_list, error, '.','MarkerSize',25, 'Color', data_color(j,:),'Displayname',sprintf("%s Data",method(j)));
            grid on;

            % plot styling
            set(gca,"TickLabelInterpreter",'latex','FontSize',12);
            legend('Location','northwest', 'Interpreter', 'Latex','FontSize',12);
            title('Global Truncation Error $\epsilon$ for Fixed-Step RK Methods vs Timestep $h$', 'Interpreter', 'Latex', 'FontSize',18);
            xlabel("Step size $h$ (-)", 'Interpreter', 'Latex','FontSize',16);
            ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
            fontsize(scale=1.5)
            axis([1e-5 10 2e-13 30])
        end 
    end
end