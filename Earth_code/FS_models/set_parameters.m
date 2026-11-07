function [a, r_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile] = set_parameters()

a = 6370000;  %planetary radius, in [m]
r_core = 3480000;  %border of core
r_mantle = 6290000;  %border of mantle

% %HOMOGENEOUS DATA [KELVIN SPHERE]
% rho_profile = [5500, 5500, 5500];           %in [kg / m^3]
% mu_profile = [145e9, 145e9, 145e9];   %in [Pa]
% lambda_profile = [352e9, 352e9, 352e9];  %in [Pa]

%PREM-LIKE DATA
rho_profile = [11000, 4800, 3200];      %in [kg / m^3]
mu_profile = [0.25e9, 190e9, 57e9];   %in [Pa]
lambda_profile = [960e9, 300e9, 70e9];  %in [Pa]

end
