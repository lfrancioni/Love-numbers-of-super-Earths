%-------------------------------------------------------------------------%
% --- Set the radii of each layer and the internal differentiation -------%
% --- profiles of the planet ---------------------------------------------%
%-------------------------------------------------------------------------%

function [a, r_in_core, r_out_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile] = set_parameters() 

%SUPER-EARTH 2M DATA
a = 7690000;          %radius, in [m]
r_in_core = 1445000;  %border of inner core, in [m]
r_out_core = 4120000;  %border of outer core, in [m]
r_mantle = 7580000;   %border of mantle, in [m]

rho_profile = [16000, 14000, 6000, 3000];      %in [kg / m^3]
mu_profile = [374e9, 0, 382e9, 57e9];   %in [Pa]
lambda_profile = [2200e9, 2085e9, 531e9, 70e9];  %in [Pa]

end
