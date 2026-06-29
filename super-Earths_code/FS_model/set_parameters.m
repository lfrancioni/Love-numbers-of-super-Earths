%-------------------------------------------------------------------------%
% --- Set the radii of each layer and the internal differentiation -------%
% --- profiles of the planet ---------------------------------------------%
%-------------------------------------------------------------------------%

function [a, r_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile] = set_parameters()

a = 11120000;  %planetary radius, in [m]

r_core = 5820000;  %border of core, in [m]
r_mantle = 11000000;  %border of mantle, in [m]

% %KELVIN SPHERE DATA
% rho_profile = [5500, 5500, 5500];           %in [kg / m^3]
% mu_profile = [1.45e11, 1.45e11, 1.45e11];   %in [Pa]
% lambda_profile = [3.5e11, 3.5e11, 3.5e11];  %in [Pa]

% %PREM-LIKE DATA
% rho_profile = [11000, 4800, 3200];      %in [kg / m^3]
% mu_profile = [2.5e8, 1.9e11, 5.7e10];   %in [Pa]
% lambda_profile = [9.6e11, 3e11, 7e10];  %in [Pa]

%SUPER-EARTH 8M DATA
rho_profile = [21000, 8000, 4500];      %in [kg / m^3]
mu_profile = [4.72e11, 7.70e11, 5.7e10];   %in [Pa]
lambda_profile = [76e11, 11.8e11, 7e10];  %in [Pa]


end