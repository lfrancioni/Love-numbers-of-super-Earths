function f_step = system_fluid(r_step, y_step, n, G, rho, rho_profile, ...
    r_in_core, r_out_core, r_mantle)

%elastic parameters in the current step
g = get_g(r_step, r_in_core, r_out_core, r_mantle, rho_profile, G);

%matrix 2x2
M_f = [(4*pi*G*rho)/g - (n+1)/r_step, 1; 
    ...
    8*pi*G*rho*(n-1)/(r_step*g), (n-1)/r_step - (4*pi*G*rho)/g];

f_step = M_f * (y_step)';

end