% Function to plot an orbit, comparing an analytical solution and three
% Runge-Kutta methods.
function plot_embedded_orbit(growth_limit)
    arguments
        growth_limit (1,:) double = 1.5   %Anything from 1.5 - 10 
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
    % set the Butcher tableaus
    BT_struct = init_BT("Darmond Prince")
    % set the reference timestep size, initial value, and time range
    h_ref = 1e-4;
    X0 = V_list(1,:);
    tspan = [0 35];
    % create a wrapper function
    rate_func_in = @(t,V) gravity_rate_func(t,V,orbit_params);

    % compute the RK methods
    [t_4o,fourth_order] = explicit_RK_variable_step_integration(rate_func_in,tspan, ...
        X0,h_ref,BT_struct);

    % plot the orbit path around the sun
    figure(1)
    axis equal; axis square;
    axis([-4, 10, -8, 6]);
    hold on;
    grid on;
    plot(0,0,'ro','markerfacecolor',[1 0.65 0],'color',[1 0.65 0],'markersize',12,DisplayName="Sun");
    plot(V_list(:,1),V_list(:,2),'k', LineWidth=4, DisplayName="Analytical Solution");
    plot(fourth_order(:,1),fourth_order(:,2),'-.',LineWidth=6,Color=[0.5 1 0.5],DisplayName='Fourth Order');
    plot(euler(:,1),euler(:,2),'--',LineWidth=5,Color=[1 0.5 0.5],DisplayName='Forward Euler');
    plot(midpoint(:,1),midpoint(:,2),':',LineWidth=4,Color=[0.5 0.5 1],DisplayName='Explicit Midpoint');
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
    plot(t_4o,fourth_order(:,1),'-.',LineWidth=4,Color=[0.5 1 0.5]);
    plot(t_euler,euler(:,1),'--',LineWidth=4,Color=[1 0.5 0.5]);
    plot(t_midpoint,midpoint(:,1),':',LineWidth=3,Color=[0.5 0.5 1]);
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
    plot(t_4o,fourth_order(:,2),'-.',LineWidth=4,Color=[0.5 1 0.5],DisplayName='Fourth Order');
    plot(0,0,'.w',DisplayName="")
    plot(t_euler,euler(:,2),'--',LineWidth=4,Color=[1 0.5 0.5],DisplayName='Forward Euler');
    plot(0,0,'.w',DisplayName="")
    plot(t_midpoint,midpoint(:,2),':',LineWidth=3,Color=[0.5 0.5 1],DisplayName='Explicit Midpoint');
    set(gca,"TickLabelInterpreter",'latex','FontSize',14);
    title("$Y$-values vs Time", "Interpreter","Latex",'FontSize',18);
    xlabel('Time (-)', 'Interpreter', 'Latex');
    ylabel('$Y$ Position (-)', 'Interpreter', 'Latex');

    my_legend = legend();
    my_legend.Layout.Tile = 'south';
    my_legend.ItemTokenSize(1) = 50;
    set(my_legend,"interpreter","latex","Fontsize",16,'Orientation', 'Horizontal');
end