%-------------------------------------------------------------------------%
% --- Set the numerical integration parameters ---------------------------%
%-------------------------------------------------------------------------%
% --- INTEGRATE FROM NEAR THE CENTER AND NOT FROM CMB DUE TO -------------%
% --- INSTABILITIES AND MISSING DEGREES LNs ----------------,,------------%
%-------------------------------------------------------------------------%


function [step, r0, r, N] = set_integration(a)

step = 10000;    %integration step, in [m]
r0 = step;       %integration starting point
r = r0:step:a;   %radial array of steps
N = length(r);    %n. of steps

end