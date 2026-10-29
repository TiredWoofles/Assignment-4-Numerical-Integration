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
function error = plot_global_truncation(run, error, var_test, func, sol, t_span, num_iter)
    arguments
        run (1,:) double = 1
        error (1,:) double = 0
        var_test (1,:) string = "Fun Count"
        func (1,:) function_handle = @gravity_wrapper
        sol (1,:) function_handle = @compute_planetary_motion
        t_span (1,:) double = [0, 30]
        num_iter (1,:) double = 50
    end

    % initialize orbit parameters
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    % set methods
    method = ["Fehlberg"];
    data_color = [0.65 .25 0; 0.25 0 0.65; 0.25 0.65 0];
    fit_color = [1 .65 0.35; 0.65 0.35 1; 0.75 1 0];

    % set initial position and velocity
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    % create wrapper functions for the gravity function and the analytical
    % solution
    sol = @(t) compute_planetary_motion(t, V0, orbit_params);

    % set initial x guess, timestep size, a filter based on the timestep
    % size (for calculating p), and loop iterations
    X0 = sol(t_span(1));
    h_list = logspace(-5,1,num_iter);
    h_filter = h_list <= 1 & h_list >= 1e-3;
    log_h_list = log10(h_list(h_filter))';
    iter = 0

    % set up the figure outside of the for-loop & if-statements
    figure(1)

    % Loop through the three RK methods
    for j = 1:length(method)
        if run
            % set up the Butcher tableau corresponding to the current method
            % loop, as well as empty lists for storing data this loop
            BT_struct = init_BT(method(j));
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
        end
    
        % set error as a log spread
        log_error = log10(error(h_filter));
    
        % Calculate p-value, plot a fit line, and plot data for a
        % comparison of error and function count.
        if var_test == "Fun Count"
            log_evals = log10(fun_count(h_filter));
            slope = polyfit(log_evals,log_error, 1);
            p = slope(1)
    
            fit_line = 10.^polyval(slope, log10(h_list));
            loglog(h_list, fit_line, '-', 'Linewidth', 7, 'Color', fit_color(j,:), 'Displayname',sprintf("Fixed Step %s Fit",method));
            hold on;
            grid on;
            loglog(fun_count, error, '.', 'MarkerSize',25, 'Color', data_color(j,:),'Displayname',sprintf("Fixed Step %s Data",method));
            
            % plot styling
            title('Global Truncation Error $\epsilon$ for RK  vs Efficiency', 'Interpreter', 'Latex', 'FontSize',18);
            legend('Location','northeast', 'Interpreter', 'Latex','FontSize',12);
            xlabel("Number of Rate Function Calls (-)", 'Interpreter', 'Latex','FontSize',16);
            ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
            set(gca,"TickLabelInterpreter",'latex');
            fontsize(scale=1.5)
            axis([3 4e6 1e-10 1e5])
    
        % Alternatively, calculate p-value, plot a fit line, and plot data 
        % for a comparison of error and step size.
        elseif var_test == "Step Size"
            slope = polyfit(log_h_list, log_error, 1);
            fit_line = 10.^polyval(slope, log10(h_list));
            p = slope(1);
    
            loglog(h_list, fit_line, '-', 'Linewidth', 7, 'Color', fit_color(j,:), 'Displayname',sprintf("Fixed Step %s Fit",method));
            hold on; % call hold on after calling loglog
            loglog(h_list, error, '.','MarkerSize',20, 'Color', data_color(j,:),'Displayname',sprintf("Fixed Step %s Data",method(j)));
            grid on;

            % plot styling
            title('Global Truncation Error $\epsilon$ for Three RK Methods vs Timestep $h$', 'Interpreter', 'Latex', 'FontSize',18);
            legend('Location','northwest', 'Interpreter', 'Latex','FontSize',12);
            xlabel("Step size $h$ (-)", 'Interpreter', 'Latex','FontSize',16);
            ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
            set(gca,"TickLabelInterpreter",'latex','FontSize',12);
            fontsize(scale=2)
            axis([h_list(1) h_list(end) 1e-12 1e3])
        end 
    end
end