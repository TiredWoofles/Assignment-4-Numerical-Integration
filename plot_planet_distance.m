% Function to plot distance from the sun as a function of step size for the
% orbit mapping function.
%INPUTS:
% method: the method used by the variable-step integrator
%OUTPUTS:
% None
function plot_planet_distance(method)
    % set optional inputs
    arguments
        method (1,:) string = "Fehlberg"
    end

    % set initial parameters
    func = @gravity_wrapper;
    sol = @compute_planet_wrapper;
    t_span = [0, 1000];
    XA = sol(t_span(end))';
    h_ref = 3;
    [BT_struct,p] = init_BT(method);
    error_desired = 0.001;
    growth_limit = 1.5;
    sun = [0, 0];

    % loop through different step sizes and calculate the truncation
    % errors for each one
    [t_list,x_list] = explicit_RK_variable_step_integration ...
        (func,t_span,XA,h_ref,BT_struct,p,error_desired,growth_limit);
    
    % store location coordinates from x & timesteps (by calculating from
    % t_list)
    location = x_list(:,1:2);
    h_list = diff(t_list)';

    % loop through location values and individually use norm to calculate
    % the distance from the sun
    % there's probably a better way to do this
    iter = length(x_list);
    distance = zeros(iter,1);
    for i = 1:iter
        distance(i) = norm(location(i,:)-sun);
    end
    distance = distance(2:end);

    % plot distance vs step size
    clf
    figure(1)
    scatter(h_list,distance,'.',LineWidth=2)
    grid on;
    
    % plot styling
    set(gca,"TickLabelInterpreter",'latex','FontSize',12);
    title('Distance from Sun, by Step Size', 'Interpreter', 'Latex', 'FontSize',18);
    % legend('Location','northeast', 'Interpreter', 'Latex','FontSize',12);
    xlabel("Step Size $h$ (-)", 'Interpreter', 'Latex','FontSize',16);
    ylabel("Distance (-)", 'Interpreter', 'Latex','FontSize',16);
    fontsize(scale=1.5)
end