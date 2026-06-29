function [a, r_in_core, r_out_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile] = set_parameters() 

a = 6370000;          %radius, in [m]

r_in_core = 1220000;
r_out_core = 3480000;
r_mantle = 6340000;

%PREM-LIKE DATA (see approach from: Poulsen 2009)
rho_profile = [12900, 11000, 4800, 3200];        %in [kg / m^3]
mu_profile = [1.7e11, 0, 1.9e11, 5.7e10];        %in [Pa]
lambda_profile = [12.7e11, 9.4e11, 3e11, 7e10];  %in [Pa]

end