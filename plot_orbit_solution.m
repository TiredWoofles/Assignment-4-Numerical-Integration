% Function to plot an orbit, comparing an analytical solution and three
% Runge-Kutta methods.
%INPUTS (all optional):
% method: An array of strings corresponding to the Butcher tableaus for the
%         methods desired, in order.
%OUTPUTS:
% None
function plot_orbit_solution(method)
    % set optional inputs
    arguments
        method (1,:) string = ["Fourth Order", "Forward Euler", "Explicit Midpoint"]
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
    color = [0.5 1 0.5; 1 0.5 0.5; 0.5 0.5 1];
    weight = [5, 4, 3];
    style = ["-.","--",":"];
    % set the reference timestep size, initial value, and time range
    h_ref = 1e-4;
    X0 = V_list(1,:);
    tspan = [0 35];
    % create a wrapper function
    rate_func_in = @(t,V) gravity_rate_func(t,V,orbit_params);

    for i = 1:length(method)
        % on the first loop only, plot the analytical solution first and
        % set up plot styling
        if i == 1
            clf
            % plot the orbit onto figure 1
            figure(1)
            axis equal; axis square;
            axis([-4, 10, -8, 6]);
            hold on;
            grid on;
            plot(0,0,'ro','markerfacecolor',[1 0.65 0],'color',[1 0.65 0],'markersize',12,DisplayName="Sun");
            plot(V_list(:,1),V_list(:,2),'k', LineWidth=4, DisplayName="Analytical Solution");
            set(gca,"TickLabelInterpreter",'latex','FontSize',18);
            title("Three Runge-Kutta Methods vs Analytical Solution", "Interpreter","Latex",'FontSize',24);
            xlabel('$X$ Position (-)', 'Interpreter', 'Latex','FontSize',20);
            ylabel('$Y$ Position (-)', 'Interpreter', 'Latex','FontSize',20);
                 
            % plot x-y value comparisons onto two figure 2 tiles
            figure(2)
            t = tiledlayout(2, 1, "tilespacing", "compact", "padding", "compact");
            title(t,"Three Runge-Kutta Methods vs Analytical Solution","Interpreter","Latex",'FontSize',34);
            
            nexttile(1)
            hold on;
            grid on;
            plot(t_range,V_list(:,1),'k', LineWidth=3);
            set(gca,"TickLabelInterpreter",'latex','FontSize',14);
            title("$X$-values vs Time", "Interpreter","Latex",'FontSize',18);
            xlabel('Time (-)', 'Interpreter', 'Latex');
            ylabel('$X$ Position (-)', 'Interpreter', 'Latex');
            
            nexttile(2)
            hold on;
            grid on;
            plot(t_range,V_list(:,2),'k', LineWidth=3, DisplayName="Analytical Solution");
            plot(0,0,'.w',DisplayName="")
            set(gca,"TickLabelInterpreter",'latex','FontSize',14);
            title("$Y$-values vs Time", "Interpreter","Latex",'FontSize',18);
            xlabel('Time (-)', 'Interpreter', 'Latex');
            ylabel('$Y$ Position (-)', 'Interpreter', 'Latex');
        end

        BT_struct = init_BT(method(i));
        % compute the RK methods
        [t_list,x_list] = explicit_RK_fixed_step_integration(rate_func_in,tspan, ...
            X0,h_ref,BT_struct);

        % plot the orbit path around the sun
        figure(1)
        plot(x_list(:,1),x_list(:,2),style(i),LineWidth=weight(i),Color=color(i,:),DisplayName=sprintf("%s",method(i)));
        legend("location","eastoutside","interpreter","latex","Fontsize",16)

        % plot the x and y values for each time step
        figure(2);
        % x-value plot:
        nexttile(1);
        plot(t_list,x_list(:,1),style(i),LineWidth=weight(i),Color=color(i,:));
        % y-value plot:
        nexttile(2)
        plot(t_list,x_list(:,2),style(i),LineWidth=weight(i),Color=color(i,:),DisplayName=method(i));
        plot(0,0,'.w',DisplayName="")
        my_legend = legend();
        my_legend.Layout.Tile = 'south';
        my_legend.ItemTokenSize(1) = 50;
        set(my_legend,"interpreter","latex","Fontsize",16,'Orientation', 'Horizontal');
    end
end