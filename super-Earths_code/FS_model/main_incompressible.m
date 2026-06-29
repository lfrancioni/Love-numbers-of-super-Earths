%-------------------------------------------------------------------------%
% --- Evaluate LNs of harmonic degree n and 'Love_type' type -------------%
% --- for incompressible planet, in a range of the specified -------------%
% --- elastic parameter --------------------------------------------------%
%-------------------------------------------------------------------------%

function [H, K, L] = main_incompressible(n, Love_type, mu_mantle)

%PARAMETERS SETTINGS-------------------------------------------------------
G = 6.67 * 10^(-11);  %in [N m^2 / kg^2]

[a, r_core, r_mantle, rho_profile, mu_profile, lambda_profile] = ...
    set_parameters();
%lambda_profile will be unused in the incompressible case

[step, r0, r, N] = set_integration(a, r_core);

%select the elastic parameter to test in a range
mu_profile(2) = mu_mantle;


%VARIABLES SCALING --------------------------------------------------------
[a, G, r_core, r_mantle, rho_profile, mu_profile, lambda_profile, ...
    step, r0, r] = normalize(a, G, r_core, r_mantle, rho_profile, ...
    mu_profile, lambda_profile, step, r0, r);


%STARTING SOLUTIONS IN r0--------------------------------------------------
rho_0 = get_parameter(r0, r_core, r_mantle, rho_profile);
mu_0 = get_parameter(r0, r_core, r_mantle, mu_profile);
g_0 = get_g(r0, r_core, r_mantle, rho_profile, G);

%Starting solution 1
y1_s0_1 = n*(r0^(n+1))/(2*(2*n + 3));
y2_s0_1 = (n+3)*(r0^(n+1))/(2*(n+1)*(2*n+3));
y3_s0_1 = (n*rho_0*g_0*r0 + 2*mu_0*(n^2 -n -3))*(r0^n)/(2*(2*n +3));
y4_s0_1 = n*(n+2)*mu_0*(r0^n)/((2*n +1)*(2*n +3));
y5_s0_1 = 0;
y6_s0_1 = 2*pi*G*rho_0*n*(r0^(n+1))/(2*n +3);

%Starting solution 2
y1_s0_2 = r0^(n-1);
y2_s0_2 = (r0^(n-1))/n;
y3_s0_2 = (rho_0*g_0*r0 + 2*(n-1)*mu_0)*(r0^(n-2));
y4_s0_2 = 2*(n-1)*mu_0*(r0^(n-2))/n;
y5_s0_2 = 0;
y6_s0_2 = 4*pi*G*rho_0*(r0^(n-1));

%Starting solution 3
y1_s0_3 = 0;
y2_s0_3 = 0;
y3_s0_3 = rho_0*(r0^n);
y4_s0_3 = 0;
y5_s0_3 = r0^n;
y6_s0_3 = (2*n +1)*(r0^(n-1));

%column arrays of starting sets
Y_s0_1 = [y1_s0_1; y2_s0_1; y3_s0_1; y4_s0_1; y5_s0_1; y6_s0_1];
Y_s0_2 = [y1_s0_2; y2_s0_2; y3_s0_2; y4_s0_2; y5_s0_2; y6_s0_2];
Y_s0_3 = [y1_s0_3; y2_s0_3; y3_s0_3; y4_s0_3; y5_s0_3; y6_s0_3];

%matrix of starting sets
Y_s0 = [Y_s0_1, Y_s0_2, Y_s0_3];


%NUMERICAL INTEGRATION VIA RUNGE KUTTA 4-----------------------------------
matrix4 = [0, 0, 0, 0, 0; 1/2, 1/2, 0, 0, 0; 1/2, 0, 1/2, 0, 0; ...
    1, 0, 0, 1, 0; 0, 1/6, 1/3, 1/3, 1/6];

Y_final = RK4_solid_incom(2, N, step, r, n, G, rho_profile, ...
    mu_profile, Y_s0, matrix4, r_core, r_mantle);


%SURFACE BOUNDARY CONDITIONS TO PRODUCE THE LINEAR COMBINATION OF y
%SOLUTION -----------------------------------------------------------------
g_N = get_g(a, r_core, r_mantle, rho_profile, G);

My = [Y_final(3,1), Y_final(3,2), Y_final(3,3); ...
    Y_final(4,1), Y_final(4,2), Y_final(4,3); ...
    Y_final(6,1), Y_final(6,2), Y_final(6,3)];

switch Love_type
    case 'Tidal'
        CC = [0; 0; (2*n + 1)/a]; 
    case 'Load'
        CC = [((2*n + 1)/3)*(3*g_N/(4*pi*a*G)); 0; (2*n +1)/a];
    case 'Shear'
        CC = [0; (2*n + 1)*g_N/(4*pi*G*a*n*(n+1)); 0];                                              
end

Const = My \ CC;  %Const = (C1, C2, C3) for the correct linear combination 
                  %at the surface


%LOVE NUMBERS--------------------------------------------------------------
y1 = Const(1)*Y_final(1,1) + Const(2)*Y_final(1,2) + Const(3)*Y_final(1,3);
y2 = Const(1)*Y_final(2,1) + Const(2)*Y_final(2,2) + Const(3)*Y_final(2,3);
y5 = Const(1)*Y_final(5,1) + Const(2)*Y_final(5,2) + Const(3)*Y_final(5,3);

switch Love_type
    case 'Shear'     %no potential --> signs ok
        H = y1*g_N;
        L = y2*g_N;
        K = y5;
    otherwise        %potential with sign --> signs to be changed
        H = -y1*g_N;
        L = -y2*g_N;
        K = y5 - 1;
end

% %Analytical LN for generic degree n, incompressible homogeneous case
% H_analitic = (2*n + 1)/(2*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));
% K_analitic = 3/(2*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));
% L_analitic = 3/(2*n*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));

end