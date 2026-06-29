%-------------------------------------------------------------------------%
% --- Set the radii of each layer and the internal differentiation -------%
% --- profiles of the planet ---------------------------------------------%
%-------------------------------------------------------------------------%

function [a, r_in_core, r_out_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile] = set_parameters() 

a = 7690000;          %radius, in [m]

r_in_core = 1445000;  %border of inner core, in [m]
r_out_core = 4120000;  %border of outer core, in [m]
r_mantle = 7580000;   %border of mantle, in [m]

% %PREM-LIKE DATA
% rho_profile = [12900, 11000, 4800, 3200];        %in [kg / m^3]
% mu_profile = [1.7e11, 0, 1.9e11, 5.7e10];        %in [Pa]
% lambda_profile = [12.7e11, 9.4e11, 3e11, 7e10];  %in [Pa]

%SUPER-EARTH 2M DATA
rho_profile = [16000, 14000, 6000, 3000];      %in [kg / m^3]
mu_profile = [3e12, 0, 3.82e11, 5.7e10];   %in [Pa]
lambda_profile = [22e11, 20.85e11, 5.31e11, 7e10];  %in [Pa]

end