function [step, r0, r, N] = set_integration(a, r_core)

step = 10000;        %integration step, in [m]
r0 = r_core - step;  %starting point for integration
r = r0:step:a;       %array of radial steps
N = length(r);       %n. of steps

end