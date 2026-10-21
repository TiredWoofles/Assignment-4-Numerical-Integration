function plot_global_truncation(step_func, var_test, rate_func, sol, t_span, num_iter)
    arguments
        step_func (1,:) function_handle
        var_test (1,:) string = "Step Size"
        rate_func (1,:) function_handle = @gravity_rate_func
        sol (1,:) function_handle = @compute_planetary_motion
        t_span (1,:) double = [0, 2*pi]
        num_iter (1,:) double = 100
    end

    X0 = sol(t_span(1));
    h_list = logspace(-5,1,num_iter);
    h_filter = h_list <= 1e-2;
    log_h_list = log10(h_list(h_filter))';
    iter = 0;
    
    figure(1)

    % Calculation
    error = zeros(num_iter,1);
    fun_count = zeros(num_iter,1);
    for i = 1:num_iter
        h = h_list(i);
        [~, x_list,~, iter] = step_func(rate_func, t_span, X0, h);  
        x_true = sol(t_span(end));              % calculate error
        error(i) = norm(x_list(end,:).' - x_true);
        fun_count(i) = iter;                    % save iterationsx
        iter = iter + 1
    end

    log_error = log10(error(h_filter));

    if var_test == "Fun Count"
        %%%%%% Fun Count P-values
        log_evals = log10(fun_count(h_filter));
        slope = polyfit(log_evals,log_error, 1);
        fit_line = 10.^polyval(slope, log10(h_list));
        p = implicit_slope(1)

        loglog(h_list, fit_line, '-', 'Linewidth', 7, 'Color', [0.5 1 0.5], 'Displayname',"Exp. Mid. Fit");
        loglog(fun_count, error, '.', 'MarkerSize',25, 'Color', [0.5 1 0.5],'Displayname',"Implicit Midpoint");
        
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
        slope = polyfit(log_h_list, log_implicit_error, 1);
        fit_line = 10.^polyval(slope, log10(h_list));
        p = implicit_slope(1)

        loglog(h_list, fit_line, '-', 'Linewidth', 7, 'Color', [0.5 1 0.5], 'Displayname',"Exp. Mid. Fit");
        hold on; % call hold on after calling loglog
        loglog(h_list, error, '.','MarkerSize',25, 'Color', [0.5 1 0.5],'Displayname',"Implicit Midpoint");
    
        % plot styling
        title('Global Truncation Error $\epsilon$ for Three RK Methods vs Timestep $h$', 'Interpreter', 'Latex', 'FontSize',18);
        legend('Location','northwest', 'Interpreter', 'Latex','FontSize',12);
        xlabel("Step size $h$ (-)", 'Interpreter', 'Latex','FontSize',16);
        ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
        set(gca,"TickLabelInterpreter",'latex');
        fontsize(scale=1.5)
        axis([h_list(1) h_list(end) 1e-10 1e5])
    end 
    
end