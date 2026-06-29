%-------------------------------------------------------------------------%
% --- Normalize all the variables of the code ----------------------------%
%-------------------------------------------------------------------------%

function [a, G, r_core, r_mantle, rho_profile, mu_profile, lambda_profile, ...
    passo, r0, r] = normalize(a, G, r_core, r_mantle, rho_profile, ...
    mu_profile, lambda_profile, passo, r0, r)

%Scaling parameters
%for lenghts
L_scale = a;

%for densities
Rho_scale = (rho_profile(1) + rho_profile(2) + rho_profile(3)) / 3;

%for times
T_scale = 1 / sqrt(Rho_scale * G * pi);


%Scaling
a = a/L_scale;
G = 1/pi;
r_core = r_core/L_scale;
r_mantle = r_mantle/L_scale;

rho_profile = rho_profile / Rho_scale; %in [kg / m^3]
mu_profile = (T_scale^2 / (L_scale^2 * Rho_scale)) * mu_profile; %in [Pa]
lambda_profile = (T_scale^2 / (L_scale^2 * Rho_scale)) * lambda_profile; %in [Pa]

passo = passo/L_scale;
r0 = r0/L_scale;
r = r / L_scale;

end