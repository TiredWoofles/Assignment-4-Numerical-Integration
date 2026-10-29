% Function to plot global truncation error against step size or function
% count for various Runge-Kutta method.
%INPUTS (all optional):
% var_test: a string that calls function count or step size calculations;
%           defaults to step size
% func: a function handle that defines a rate function to approximate
% sol: a function handle that defines the solution of func
% t_span: a two-value list containing the first and final values desired of
%         the relevant time span
% num_iter: the number of iterations/datapoints desired
% h_ref: a starting h-value (timestep size)
% growth_limit: the maximum factor that h is allowed to grow by
%OUTPUTS:
% None
function plot_embedded_global_truncation(method, var_test, func, sol, ...
    t_span, error_span, num_iter, h_ref, growth_limit)
    arguments
        method (1,:) string = ["Fehlberg"]
        var_test (1,:) string = "Step Size"
        func (1,:) function_handle = @gravity_wrapper
        sol (1,:) function_handle = @compute_planet_wrapper
        t_span (1,:) double = [0, 30]
        error_span (1,:) double = [-13, 1]
        num_iter (1,:) double = 50
        h_ref (1,:) double = 1e-2;
        growth_limit = 8;
    end

    % set method
    data_color = [0.5 0 0; 0.25 0.25 0.65; 0.25 0.65 0];
    fit_color = [1 0.45 0.45; 0.65 0.65 1; 0.6 1 0.6];

    % set initial x guess, timestep size, a filter based on the timestep
    % size (for calculating p), and loop iterations
    X0 = sol(t_span(1))';
    e_list = logspace(error_span(1),error_span(2),num_iter);
    e_list = e_list(end:-1:1);
    iter = 0;

    % Loop through the three RK method
    for j = 1:length(method)
        % set up the Butcher tableau corresponding to the current method
        % loop: BT, p-value, and filter bounds
        [BT_struct,p,bounds] = init_BT(method(j),var_test);
        % set up empty lists for storing data this loop
        error = zeros(num_iter,1);
        fun_count = zeros(num_iter,1);
        failure_list = zeros(num_iter,1);
        h_avg_list = zeros(num_iter,1);

        % loop through different step sizes and calculate the truncation
        % errors for each one
        for i = 1:num_iter
            error_desired = e_list(i);
            [~,x_list,h_avg, evals,failure_rate] = explicit_RK_variable_step_integration ...
                (func,t_span,X0,h_ref,BT_struct,p,error_desired,growth_limit);
            x_true = sol(t_span(end));
            
            error(i) = norm(x_list(end,:).' - x_true);
            h_avg_list(i) = h_avg;
            fun_count(i) = evals;
            failure_list(i) = failure_rate;
            
            % % print iterations for troubleshooting
            iter = iter + 1;
        end
        

        % set up the figure
        figure(1)
    
        % Calculate p-value, plot a fit line, and plot data for a
        % comparison of error and function count.
        if var_test == "Fun Count"
            % set up a filter for the linear portion of the method
            filter = fun_count >= bounds(1) & fun_count <= bounds(2);
            log_error = log10(error(filter));
            log_evals = log10(fun_count(filter));

            % calculate the p-value (aka the slope) and fit line
            slope = polyfit(log_evals,log_error, 1);
            fit_line = 10.^polyval(slope, log10(fun_count));
            p = slope(1)
    
            % plot the data
            loglog(fun_count, fit_line, '-', 'Linewidth', 7, 'Color', fit_color(j,:), 'Displayname',sprintf("%s Fit", method(j)));
            hold on; % call hold on AFTER loglog
            loglog(fun_count, error, '.', 'MarkerSize',25, 'Color',data_color(j,:),'Displayname',sprintf("%s Data Points",method(j)));
            grid on;
            
            % plot styling
            set(gca,"TickLabelInterpreter",'latex','FontSize',12);
            title('Global Truncation Error $\epsilon$ for RK method vs Efficiency', 'Interpreter', 'Latex', 'FontSize',18);
            legend('Location','northeast', 'Interpreter', 'Latex','FontSize',12);
            xlabel("Number of Rate Function Calls (-)", 'Interpreter', 'Latex','FontSize',16);
            ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
            fontsize(scale=1.5)
            axis([fun_count(1) fun_count(end) 1e-12 1e3])
    
        % Alternatively, calculate p-value, plot a fit line, and plot data 
        % for a comparison of error and step size.
        elseif var_test == "Step Size"
            % set up a filter for the linear portion of the method
            filter = h_avg_list >= bounds(1) & h_avg_list <= bounds(2);
            log_error = log10(error(filter));
            log_h_list = log10(h_avg_list(filter));

            % calculate the p-value (aka the slope) and fit line
            slope = polyfit(log_h_list, log_error, 1);
            fit_line = 10.^polyval(slope, log10(h_avg_list));
            p = slope(1)

            % plot the data
            loglog(h_avg_list, fit_line, '-', 'Linewidth', 7, 'Color', fit_color(j,:), 'Displayname',sprintf("Adaptive Step %s Fit",method(j)));
            hold on; % call hold on AFTER loglog
            loglog(h_avg_list, error, '.','MarkerSize',25, 'Color', data_color(j,:),'Displayname',sprintf("Adaptive Step %s Method",method(j)));
            grid on;

            % plot styling
            set(gca,"TickLabelInterpreter",'latex','FontSize',14);
            title('Runge-Kutta Global Truncation Error $\epsilon$ vs Timestep $h$', 'Interpreter', 'Latex', 'FontSize',18);
            legend('Location','northwest', 'Interpreter', 'Latex','FontSize',8);
            xlabel("Average Step Size $h_{avg}$ (-)", 'Interpreter', 'Latex','FontSize',16);
            ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
            fontsize(scale=2)
            axis([1e-4 1e1 1e-12 1e2])
        end 
    end
end