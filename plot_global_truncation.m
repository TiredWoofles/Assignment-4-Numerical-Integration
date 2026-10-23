% Function to plot global truncation error against step size or function
% count for various Runge-Kutta methods.
%INPUTS:
% var_test: a string that calls function count or step size calculations;
%           defaults to step size
% rate_func: a function handle that defines 
%V: the vector of the position and velocity of the planet
%   V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
%   orbit_params.m_sun is the mass of the sun
%   orbit_params.m_planet is the mass of the planet
%   orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
%   dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function error = plot_global_truncation(var_test, rate_func, sol, t_span, num_iter)
    arguments
        % error (:,1) double = 0
        var_test (1,:) string = "Step Size"
        rate_func (1,:) function_handle = @gravity_rate_func
        sol (1,:) function_handle = @compute_planetary_motion
        t_span (1,:) double = [0, 2*pi]
        num_iter (1,:) double = 50
    end
    % set orbit parameters
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    % set Butcher tableaus
    BT_4o = struct();
    BT_4o.A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
    BT_4o.B = [1/6 1/3 1/3 1/6];
    BT_4o.C = [0; 0.5; 0.5; 1];
    BT_euler = struct();
    BT_euler.A = [0];
    BT_euler.B = [1];
    BT_euler.C = [0];
    BT_midpoint = struct();
    BT_midpoint.A = [0 0; 0.5 0];
    BT_midpoint.B = [0 1];
    BT_midpoint.C = [0; 0.5];

    methods = ["Fourth Order", "Forward Euler", "Explicit Midpoint"];
    fit_color = [0.65 0 0; 0 0 0.65; 0 0.65 0];
    data_color = [1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5];
    BT = [BT_4o, BT_euler, BT_midpoint];

    % set initial position and velocity
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    % create wrapper functions for the gravity function and the analytical
    % solution
    rate_func = @(t, V) gravity_rate_func(t, V, orbit_params);
    sol = @(t) compute_planetary_motion(t, V0, orbit_params);

    % set initial x guess, timestep size, and a filter based on the time
    % step size (for calculating p)
    X0 = sol(t_span(1));
    h_list = logspace(-5,1,num_iter);
    h_filter = h_list <= 1 & h_list >= 1e-3;
    log_h_list = log10(h_list(h_filter))';
    iter = 0;

    % set up the figure outside of the for loop & if-statements
    figure(1)

    % Loop through the three RK methods
    for j = 1:length(methods)
        BT_midpoint = BT(j);
        error = zeros(num_iter,1);
        fun_count = zeros(num_iter,1);

        % loop through different step sizes and calculate the truncation
        % errors for each one
        for i = 1:num_iter
            h = h_list(i);
            [~, x_list,~, fun_count] = explicit_RK_fixed_step_integration ...
                (rate_func,t_span,X0,h,BT_midpoint);  
            x_true = sol(t_span(end));              % calculate error
            error(i) = norm(x_list(end,:).' - x_true);
            fun_count(i) = iter;                    % save iterationsx
            iter = iter + 1
        end
    
        % set error as a logs spread
        log_error = log10(error(h_filter));
    
        % plot the graph depending on the value we're tracking
        if var_test == "Fun Count"
            %%%%%% Fun Count P-values
            log_evals = log10(fun_count(h_filter));
            slope = polyfit(log_evals,log_error, 1);
            fit_line = 10.^polyval(slope, log10(h_list));
            p = slope(1)
    
            loglog(h_list, fit_line, '-', 'Linewidth', 7, 'Color', [0.5 1 0.5], 'Displayname',"Fit");
            loglog(fun_count, error, '.', 'MarkerSize',25, 'Color', [0.5 1 0.5],'Displayname',"RK");
            
            % plot styling
            title('Global Truncation Error $\epsilon$ for RK Methods vs Efficiency', 'Interpreter', 'Latex', 'FontSize',18);
            legend('Location','northeast', 'Interpreter', 'Latex','FontSize',12);
            xlabel("Number of Rate Function Calls (-)", 'Interpreter', 'Latex','FontSize',16);
            ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
            set(gca,"TickLabelInterpreter",'latex');
            fontsize(scale=1.5)
            axis([3 4e6 1e-10 1e5])
    
        elseif var_test == "Step Size"
            %%%%%%% Step Size P-values
            slope = polyfit(log_h_list, log_error, 1);
            fit_line = 10.^polyval(slope, log10(h_list));
            p = slope(1)
    
            loglog(h_list, fit_line, '-', 'Linewidth', 3, 'Color', fit_color(j,:), 'Displayname',sprintf("%s Fit Line",methods(j)));
            hold on; % call hold on after calling loglog
            grid on;
            loglog(h_list, error, '.','MarkerSize',15, 'Color', data_color(j,:),'Displayname',sprintf("%s Data Points",methods(j)));
        
            % plot styling
            title('Global Truncation Error $\epsilon$ for Three RK Methods vs Timestep $h$', 'Interpreter', 'Latex', 'FontSize',18);
            legend('Location','northwest', 'Interpreter', 'Latex','FontSize',12);
            xlabel("Step size $h$ (-)", 'Interpreter', 'Latex','FontSize',16);
            ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
            set(gca,"TickLabelInterpreter",'latex','FontSize',12);
            fontsize(scale=1.5)
            axis([h_list(1) h_list(end) 1e-10 1e5])
        end 
    end
end