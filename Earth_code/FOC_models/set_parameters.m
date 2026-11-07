function [a, r_in_core, r_out_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile] = set_parameters() 

a = 6370000;          %radius, in [m]

r_in_core = 1220000;
r_out_core = 3480000;
r_mantle = 6340000;

%PREM-LIKE DATA (see approach from: Poulsen 2009)
rho_profile = [12900, 11000, 4800, 3200];        %in [kg / m^3]
mu_profile = [170e9, 0, 190e9, 57e9];        %in [Pa]
lambda_profile = [1270e9, 940e9, 300e9, 70e9];  %in [Pa]

end
