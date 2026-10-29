% Wrapper function for gravity_rate_func.
%INPUTS:
% t: the current timestep
% V: the initial position
%OUTPUTS:
% dVdt: the new position of the planet
function V_list = compute_planet_wrapper(t)
    % initialize orbit parameters
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    % set initial position and velocity
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];

    % compute the analytical solution
    V_list = compute_planetary_motion(t, V0, orbit_params);
end