% Wrapper function for gravity_rate_func.
%INPUTS:
% t: the current timestep
% V: the current position
%OUTPUTS:
% dVdt: the new position of the planet
function dVdt = gravity_wrapper(t, V)
    % initialize orbit parameters
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    % compute the gravity rate
    dVdt = gravity_rate_func(t, V, orbit_params);
end