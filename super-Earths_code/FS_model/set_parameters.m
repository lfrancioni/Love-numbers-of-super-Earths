%-------------------------------------------------------------------------%
% --- Set the radii of each layer and the internal differentiation -------%
% --- profiles of the planet ---------------------------------------------%
%-------------------------------------------------------------------------%

function [a, r_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile] = set_parameters()

% %SUPER-EARTH 2M DATA
%a = 7670000;  %planetary radius, in [m]
%r_core = 4130000;  %border of core, in [m]
%r_mantle = 7580000;  %border of mantle, in [m]
%
%rho_profile = [14000, 6000, 3000];      %in [kg / m^3]
%mu_profile = [374e9, 382e9, 57e9];   %in [Pa]
%lambda_profile = [2085e9, 531e9, 70e9];  %in [Pa]

% %SUPER-EARTH 4M DATA
%a = 9240000;  %planetary radius, in [m]
%r_core = 4900000;  %border of core, in [m]
%r_mantle = 9120000;  %border of mantle, in [m]
%
%rho_profile = [17000, 7000, 4000];      %in [kg / m^3]
%mu_profile = [419e9, 538e9, 57e9];   %in [Pa]
%lambda_profile = [3773e9, 789e9, 70e9];  %in [Pa]

%SUPER-EARTH 8M DATA
a = 11120000;  %planetary radius, in [m]
r_core = 5820000;  %border of core, in [m]
r_mantle = 11000000;  %border of mantle, in [m]

rho_profile = [21000, 8000, 4500];      %in [kg / m^3]
mu_profile = [472e9, 770e9, 57e9];   %in [Pa]
lambda_profile = [7604e9, 1179e9, 70e9];  %in [Pa]


end
