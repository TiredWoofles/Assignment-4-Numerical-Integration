% Function to plot an orbit using an embedded RK method
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
function plot_embedded_orbit(growth_limit, method, p)
    % set optional inputs
    arguments
        growth_limit (1,:) double = 1.5   % any value from about 1.5 to 10
        method (1,:) string = "Heun Euler"
        p (1,:) double = 2 % corresponds to the lower order of an embedded method
    end

    % define paramters for the analytical solution:
    % set the orbit mechanics
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    % set the initial location, velocity, and time range
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,35,100);
    % compute the analytical solution
    V_list = compute_planetary_motion(t_range,V0,orbit_params);

    % define parameters for the RK methods:
    % wrapper function, time range, initial guess, reference timestep,
    % Butcher tableau, error scaling, and desired error
    func = @(t,V) gravity_rate_func(t,V,orbit_params);
    tspan = [0 35];
    X0 = V_list(1,:);
    h_ref = 1e-2;
    BT_struct = init_BT(method);
    p = 2;
    error_desired = 0.3;
    
    % compute the RK methods
    [step_failure_rate, t_list,x_list,h_avg,fun_count] = explicit_RK_variable_step_integration( ...
        func,tspan,X0,h_ref,BT_struct,p,error_desired,growth_limit);

    % plot the orbit path around the sun
    clf
    figure(1)
    axis equal; axis square;
    axis([-4, 10, -8, 6]);
    hold on;
    grid on;
    plot(0,0,'ro','markerfacecolor',[1 0.65 0],'color',[1 0.65 0],'markersize',12,DisplayName="Sun");
    plot(V_list(:,1),V_list(:,2),'k', LineWidth=4, DisplayName="Analytical Solution");
    plot(x_list(:,1),x_list(:,2),':',LineWidth=3,Color=[0.5 1 0.5],DisplayName=sprintf("%s Approximation",method));
    legend("location","eastoutside","interpreter","latex","Fontsize",16)
    title("Three Runge-Kutta Methods vs Analytical Solution", "Interpreter","Latex",'FontSize',34);
    xlabel('$X$ Position (-)', 'Interpreter', 'Latex','FontSize',32);
    ylabel('$Y$ Position (-)', 'Interpreter', 'Latex','FontSize',32);
    set(gca,"TickLabelInterpreter",'latex','FontSize',18);

    % plot the x and y values for each time step
    figure(2);
    t = tiledlayout(2, 1, "tilespacing", "compact", "padding", "compact");
    title(t,"Three Runge-Kutta Methods vs Analytical Solution","Interpreter","Latex",'FontSize',34);

    % x-value plot:
    nexttile();
    hold on;
    grid on;
    plot(t_range,V_list(:,1),'k', LineWidth=3);
    plot(t_list,x_list(:,1),'-.',LineWidth=4,Color=[0.5 1 0.5]);
    set(gca,"TickLabelInterpreter",'latex','FontSize',14);
    title("$X$-values vs Time", "Interpreter","Latex",'FontSize',18);
    xlabel('Time (-)', 'Interpreter', 'Latex');
    ylabel('$X$ Position (-)', 'Interpreter', 'Latex');

    % y-value plot:
    nexttile()
    hold on;
    grid on;
    plot(t_range,V_list(:,2),'k', LineWidth=3, DisplayName='Analytical Solution');
    plot(0,0,'.w',DisplayName="")
    plot(t_list,x_list(:,2),'-.',LineWidth=4,Color=[0.5 1 0.5],DisplayName=sprintf("%s Method",method));
    plot(0,0,'.w',DisplayName="")
    set(gca,"TickLabelInterpreter",'latex','FontSize',14);
    title("$Y$-values vs Time", "Interpreter","Latex",'FontSize',18);
    xlabel('Time (-)', 'Interpreter', 'Latex');
    ylabel('$Y$ Position (-)', 'Interpreter', 'Latex');

    my_legend = legend();
    my_legend.Layout.Tile = 'south';
    my_legend.ItemTokenSize(1) = 50;
    set(my_legend,"interpreter","latex","Fontsize",16,'Orientation', 'Horizontal');
end