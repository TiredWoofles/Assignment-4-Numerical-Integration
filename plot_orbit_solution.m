%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
function [fourth_order, euler, midpoint] = plot_orbit_solution(run, fourth_order, euler, midpoint)
    arguments
        run (1,:) double
        fourth_order (:,4) double = [0 0 0 0]
        euler (:,4) double = [0 0 0 0]
        midpoint (:,4) double = [0 0 0 0]
    end
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    rate_func_in = @(t,V) gravity_rate_func(t,V,orbit_params);
    tspan = [0 15];
    h_ref = 0.00003;
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
    
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    X0 = V_list(1,:);

    if run
        [~,fourth_order] = explicit_RK_fixed_step_integration(rate_func_in,tspan, ...
            X0,h_ref,BT_4o);
        [~,euler] = explicit_RK_fixed_step_integration(rate_func_in,tspan, ...
            X0,h_ref,BT_euler);
        [~,midpoint] = explicit_RK_fixed_step_integration(rate_func_in,tspan, ...
            X0,h_ref,BT_midpoint);
    end

    axis equal; axis square;
    axis([-8, 12, -10, 10])
    hold on
    grid on
    plot(0,0,'ro','markerfacecolor',[0.75 0.45 0],'color',[0.75 0.45 0],'markersize',5,DisplayName="Sun");
    plot(V_list(:,1),V_list(:,2),'k', LineWidth=2, DisplayName="Analytical Solution");
    plot(fourth_order(:,1),fourth_order(:,2),'-.',LineWidth=4,Color=[0.5 1 0.5],DisplayName='Fourth Order');
    plot(euler(:,1),euler(:,2),'--',LineWidth=4,Color=[1 0.5 0.5],DisplayName='Forward Euler');
    plot(midpoint(:,1),midpoint(:,2),':',LineWidth=4,Color=[0.5 0.5 1],DisplayName='Implicit Midpoint');
    legend("location","eastoutside","interpreter","latex","Fontsize",16)
    title("Three Runge-Kutta Methods vs Analytical Solution", "Interpreter","Latex",'FontSize',34)
    xlabel('$X$ Position (-)', 'Interpreter', 'Latex','FontSize',32);
    ylabel('$Y$ Position (-)', 'Interpreter', 'Latex','FontSize',32);
    set(gca,"TickLabelInterpreter",'latex','FontSize',18);
end