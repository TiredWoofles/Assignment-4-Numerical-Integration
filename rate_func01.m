%Function that outputs the slope of rate_func01 at a given 
% X value and time t.
function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end