function plot_global_truncation(step_func, var_test, rate_func, sol, t_span, num_iter)
    arguments
        step_func (1,:) string
        var_test (1,:) string = "Fun Count"
        rate_func (1,:) function_handle = @rate_func02
        sol (1,:) function_handle = @solution02
        t_span (1,:) double = [0, 2*pi]
        num_iter (1,:) double = 100
    end

    X0 = sol(t_span(1));
    h_list = logspace(-5,1,num_iter);
    h_filter = h_list <= 1e-2;
    log_h_list = log10(h_list(h_filter))';
    iter = 0;
    
    figure(1)

    % Explicit Method calculations
    if step_func == "Explicit Methods" || step_func == "All"
        explicit_midpoint_error = zeros(num_iter,1);
        explicit_evals = zeros(num_iter,1);

        forward_evals = zeros(num_iter,1);
        forward_euler_error = zeros(num_iter, 1);

          for i = 1:num_iter
            h = h_list(i);
            [~, explicit_midpoint_list, ~, explicit_iter] = explicit_midpoint_fixed_step_integration(rate_func, t_span, X0, h);
            [~, forward_euler_list, ~, forward_iter] = forward_euler_fixed_step_integration(rate_func, t_span, X0, h);
            
            % calculate error
            x_true = sol(t_span(end));
            explicit_midpoint_error(i) = norm(explicit_midpoint_list(end,:).' - x_true);
            forward_euler_error(i) = norm(forward_euler_list(end,:).' - x_true);

            % save iterations
            explicit_evals(i) = explicit_iter;
            forward_evals(i) = forward_iter;
            iter = iter+1
          end
          
        log_explicit_midpoint_error = log10(explicit_midpoint_error(h_filter));
        log_forward_euler_error = log10(forward_euler_error(h_filter));

        if var_test == "Fun Count"
            %%%%% Explicit Method Fun Count P-values
            % Explicit Midpoint fit line
            log_explicit_evals = log10(explicit_evals(h_filter));
            explicit_midpoint_slope = polyfit(log_explicit_evals, log_explicit_midpoint_error, 1);
            explicit_midpoint_fit = 10.^polyval(explicit_midpoint_slope, log10(explicit_evals));
            explicit_midpoint_p = explicit_midpoint_slope(1)

            % Forward Euler fit line
            log_forward_evals = log10(forward_evals(h_filter));
            forward_euler_slope = polyfit(log_forward_evals, log_forward_euler_error, 1);
            forward_euler_fit = 10.^polyval(forward_euler_slope, log10(forward_evals));
            forward_euler_p = forward_euler_slope(1)

            % Plot explicit/fun count fit lines
            loglog(explicit_evals, explicit_midpoint_fit, '-', 'Linewidth', 7, 'Color', [0.5 1 0.5], 'Displayname',"Exp. Mid. Fit");
            hold on;
            grid on;
            loglog(forward_evals, forward_euler_fit, '-', 'Linewidth', 7, 'Color', [1 0.5 0.5], 'Displayname',"Forw. Euler Fit");
            
            % Plot explicit/step-size points
            loglog(explicit_evals, explicit_midpoint_error, '.', 'MarkerSize',25, 'Color',[0 0.65 0], 'Displayname',"Explicit Midpoint");
            hold on;
            grid on;
            loglog(forward_evals, forward_euler_error, '.', 'MarkerSize', 25, 'Color', [0.65 0 0], 'Displayname',"Forward Euler");
                
        elseif var_test == "Step Size"
            %%%%%% Explicit Method Step Size P-values
            % % Explicit Midpoint fit line
            % explicit_midpoint_step_slope = polyfit(log_h_list, log_explicit_midpoint_error, 1);
            % explicit_midpoint_fit = 10.^polyval(explicit_midpoint_step_slope, log10(h_list));
            % explicit_midpoint_p = explicit_midpoint_step_slope(1)
            % 
            % % Forward Euler fit line
            % forward_euler_step_slope = polyfit(log_h_list, log_forward_euler_error, 1);
            % forward_euler_fit = 10.^polyval(forward_euler_step_slope, log10(h_list));
            % forward_euler_p = forward_euler_step_slope(1)
            % 
            % Plot explicit/step-size fit line
            % loglog(h_list, explicit_midpoint_fit, '-', 'Linewidth', 7, 'Color', [0.5 1 0.5], 'Displayname',"Exp. Mid. Fit");
            % hold on; % call hold on after calling loglog
            % grid on;
            % loglog(h_list, forward_euler_fit, '-', 'Linewidth', 7, 'Color', [1 0.5 0.5], 'Displayname',"Forw. Euler Fit");
            
            % Plot explicit/step-size points
            loglog(h_list, explicit_midpoint_error, '.', 'MarkerSize',25, 'Color',[0 0.65 0], 'Displayname',"Explicit Midpoint");
            hold on;
            grid on;
            loglog(h_list, forward_euler_error, '.', 'MarkerSize', 25, 'Color', [0.65 0 0], 'Displayname',"Forward Euler");
        end
   
    % Implicit Midpoint calculation
    elseif step_func == "Implicit Midpoint" || step_func == "All"
        implicit_midpoint_error = zeros(num_iter,1);
        implicit_evals = zeros(num_iter,1);
        for i = 1:num_iter
            h = h_list(i);
            [~, implicit_midpoint_list,~, implicit_iter] = fixed_step_integration(rate_func, 'Implicit Midpoint', t_span, X0, h);  
            x_true = sol(t_span(end));              % calculate error
            implicit_midpoint_error(i) = norm(implicit_midpoint_list(end,:).' - x_true);
            implicit_evals(i) = implicit_iter;      % save iterations
            iter = iter+1
        end

        log_implicit_midpoint_error = log10(implicit_midpoint_error(h_filter));

        if var_test == "Fun Count"
            %%%%%% Implicit Midpoint Fun Count P-values
            log_implicit_evals = log10(implicit_evals(h_filter));
            implicit_midpoint_slope = polyfit(log_implicit_evals,log_implicit_midpoint_error, 1);
            implicit_midpoint_p = implicit_midpoint_slope(1)
            loglog(implicit_evals, implicit_midpoint_error, '.', 'MarkerSize',15, 'Color', [0.5 1 0.5],'Displayname',"Implicit Midpoint");
    
        elseif var_test == "Step Size"
            %%%%%%% Implicit Midpoint Step Size P-values
            % implicit_midpoint_slope = polyfit(log_h_list, log_implicit_midpoint_error, 1);
            % implicit_p = implicit_midpoint_slope(1)
            loglog(h_list, implicit_midpoint_error, '.','MarkerSize',15, 'Color', [0.5 1 0.5],'Displayname',"Implicit Midpoint");
        end 
    
    % Backward Euler calculations
    elseif step_func == "Backward Euler" || step_func == "All"
        backward_euler_error = zeros(num_iter,1);
        backward_evals = zeros(num_iter,1);
        for i = 1:num_iter 
            h = h_list(i);
            [~, backward_euler_list, ~, backward_iter] = fixed_step_integration(rate_func, 'Backward Euler', t_span, X0, h);    
            x_true = sol(t_span(end));              % calculate error
            backward_euler_error(i) = norm(backward_euler_list(end,:).' - x_true);
            backward_evals(i) = backward_iter;       % save iterations
            iter = iter+1
        end

        log_backward_euler_error = log10(backward_euler_error(h_filter));
        
        if var_test == "Fun Count"
            %%%%%% Backward Euler Fun Count P-values
            log_backward_evals = log10(backward_evals(h_filter));
            backward_euler_eval_slope = polyfit(log_backward_evals, log_backward_euler_error, 1);
            backward_euler_p = backward_euler_eval_slope(1)
            loglog(backward_evals, backward_euler_error, '.', 'MarkerSize',15,'Color', [1 0.5 0.5],'Displayname',"Backward Euler");
            grid on
        elseif var_test == "Step Size"
            %%%%% Backward Euler Step Size P-values
            backward_euler_slope = polyfit(log_h_list, log_backward_euler_error, 1);
            backward_euler_p = backward_euler_slope(1)
            loglog(h_list, backward_euler_error, '.', 'MarkerSize',15,'Color', [1 0.5 0.5],'Displayname',"Backward Euler");
        end
    end

    if var_test == "Fun Count"
        % plot styling
        title('Global Truncation Error $\epsilon$ for All Methods vs Efficiency', 'Interpreter', 'Latex', 'FontSize',18);
        legend('Location','northeast', 'Interpreter', 'Latex','FontSize',12);
        xlabel("Number of Rate Function Calls (-)", 'Interpreter', 'Latex','FontSize',16);
        ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
        set(gca,"TickLabelInterpreter",'latex');
        fontsize(scale=1.5)
        axis([3 4e6 1e-10 1e5])
    
   elseif var_test == "Step Size"   
        % plot styling
        title('Global Truncation Error $\epsilon$ for Four Methods vs Timesteps $h$', 'Interpreter', 'Latex', 'FontSize',18);
        legend('Location','northwest', 'Interpreter', 'Latex','FontSize',12);
        xlabel("Step size $h$ (-)", 'Interpreter', 'Latex','FontSize',16);
        ylabel("Global Truncation Error $\epsilon$ (-)", 'Interpreter', 'Latex','FontSize',16);
        set(gca,"TickLabelInterpreter",'latex');
        fontsize(scale=1.5)
        axis([h_list(1) h_list(end) 1e-10 1e5])
   end
       

end