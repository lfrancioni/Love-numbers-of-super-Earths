%-------------------------------------------------------------------------%
% --- Solve the solid ODEs system for incompressible solid layer ---------%
% --- at r_step ----------------------------------------------------------%
%-------------------------------------------------------------------------%

function f_step = system_solid_incom(r_step, y_step, n, G, rho, rho_profile, ...
    mu, r_in_core, r_out_core, r_mantle)

%parameters in the current step
g = get_g(r_step, r_in_core, r_out_core, r_mantle, rho_profile, G);

%matrix 6x3
M = [-2/r_step, n*(n+1)/r_step, 0, 0, 0, 0;
    ...
    -1/r_step, 1/r_step, 0, 1/mu, 0, 0;
    ...
    (4/r_step)*(3*mu/r_step - rho*g), -(n*(n+1)/r_step)*(6*mu/r_step - rho*g), ...
    0, n*(n+1)/r_step, -rho*(n+1)/r_step, rho;
    ...
    -(1/r_step)*(6*mu/r_step - rho*g), -(2/(r_step^2))*(1 - 2*n*(n+1))*mu, ...
    -1/r_step, -3/r_step, rho/r_step, 0;
    ...
    -4*pi*G*rho, 0, 0, 0, -(n+1)/r_step, 1;
    ...
    -4*pi*G*rho*(n+1)/r_step, 4*pi*G*rho*n*(n+1)/r_step, 0, 0, 0, (n-1)/r_step];

f_step = M * (y_step)';

end