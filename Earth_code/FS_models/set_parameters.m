function [a, r_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile] = set_parameters()

a = 6370000;  %planetary radius, in [m]
r_core = 3480000;  %border of core
r_mantle = 6290000;  %border of mantle

% %HOMOGENEOUS DATA
% rho_profile = [5500, 5500, 5500];           %in [kg / m^3]
% mu_profile = [1.45e11, 1.45e11, 1.45e11];   %in [Pa]
% lambda_profile = [3.52e11, 3.52e11, 3.52e11];  %in [Pa]

%PREM-LIKE DATA
rho_profile = [11000, 4800, 3200];      %in [kg / m^3]
mu_profile = [2.5e8, 1.9e11, 5.7e10];   %in [Pa]
lambda_profile = [9.6e11, 3e11, 7e10];  %in [Pa]

end