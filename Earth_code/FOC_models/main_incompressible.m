function [H, K, L] = main_incompressible(n, Love_type)

%planetary and integration parameters--------------------------------------
G = 6.67 * 10^(-11);  %in [N m^2 / kg^2]

[a, r_in_core, r_out_core, r_mantle, rho_profile, mu_profile, lambda_profile] = ...
    set_parameters();
%lambda_profile will be unused in this incompressible case

[passo, r0, r, N] = set_integration(a);

index_in = find(r == r_in_core);   %indexes for num integration
index_out = find(r == r_out_core);


%VARIABLES SCALING --------------------------------------------------------
[a, G, r_in_core, r_out_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile, passo, r0, r] = normalize(a, G, r_in_core, ...
    r_out_core, r_mantle, rho_profile, mu_profile, lambda_profile, ...
    passo, r0, r);


%Analytical solutions near the starting point r0 --------------------------
rho_0 = get_parameter(r0, r_in_core, r_out_core, r_mantle, rho_profile);
mu_0 = get_parameter(r0, r_in_core, r_out_core, r_mantle, mu_profile);

g_0 = get_g(r0, r_in_core, r_out_core, r_mantle, rho_profile, G);

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

%Column arrays of analytical solutions
Y_s0_1 = [y1_s0_1; y2_s0_1; y3_s0_1; y4_s0_1; y5_s0_1; y6_s0_1];
Y_s0_2 = [y1_s0_2; y2_s0_2; y3_s0_2; y4_s0_2; y5_s0_2; y6_s0_2];
Y_s0_3 = [y1_s0_3; y2_s0_3; y3_s0_3; y4_s0_3; y5_s0_3; y6_s0_3];

%matrix of starting solutions
Y_s0 = [Y_s0_1, Y_s0_2, Y_s0_3];


%NUMERICAL INTEGRATION VIA RUNGE KUTTA 4 LAYER BY LAYER--------------------
matrix4 = [0, 0, 0, 0, 0; 1/2, 1/2, 0, 0, 0; 1/2, 0, 1/2, 0, 0; 1, 0, 0, 1, 0; ...
    0, 1/6, 1/3, 1/3, 1/6];

%Inner core integration
Y_in_core = RK4_solid_incom(2, index_in, passo, r, n, G, rho_profile, ...
    mu_profile, Y_s0, matrix4, r_in_core, r_out_core, r_mantle);


%Outer core integration
%B.C. entering a fluid layer
rho_l = rho_profile(2);
g_in_core = get_g(r_in_core, r_in_core, r_out_core, r_mantle, rho_profile, G);

Frac14 = Y_in_core(4,1)/Y_in_core(4,3);
Wn1 = (Y_in_core(1,1) - Frac14 * Y_in_core(1,3));
Wn2 = (Y_in_core(5,1) - Frac14 * Y_in_core(5,3)) / g_in_core;
Wn3 = (Y_in_core(3,1) - Frac14 * Y_in_core(3,3)) / (rho_l * g_in_core);
Frac24 = Y_in_core(4,2)/Y_in_core(4,3);
Wd1 = (Y_in_core(1,2) - Frac24 * Y_in_core(1,3));
Wd2 = (Y_in_core(5,2) - Frac24 * Y_in_core(5,3)) / g_in_core;
Wd3 = (Y_in_core(3,2) - Frac24 * Y_in_core(3,3)) / (rho_l * g_in_core);

W = (Wn1 + Wn2 - Wn3) / (Wd1 + Wd2 - Wd3);

y3_0 = (Y_in_core(3,1) - Frac14 * Y_in_core(3,3)) - W * ...
    (Y_in_core(3,2) - Frac24 * Y_in_core(3,3));

y5_0 = (Y_in_core(5,1) - Frac14 * Y_in_core(5,3)) - W * ...
    (Y_in_core(5,2) - Frac24 * Y_in_core(5,3));

y6_0 = (Y_in_core(6,1) - Frac14 * Y_in_core(6,3)) - W * ...
    (Y_in_core(6,2) - Frac24 * Y_in_core(6,3));

y7_0 = y6_0 - 4*pi*G*y3_0/g_in_core;

Y_fluid = [y5_0; y7_0];
Y_out_core = RK4_fluid(index_in+1, index_out, passo, r, n, G, ...
    rho_l, rho_profile, Y_fluid, matrix4, r_in_core, r_out_core, r_mantle);


%Mantle and crust integration
g_out_core = get_g(r_out_core, r_in_core, r_out_core, r_mantle, rho_profile, G);

%B.C. exiting a fluid layer
y1_solid1 = -Y_out_core(1,1)/g_out_core;
y2_solid1 = 0;
y3_solid1 = 0;
y4_solid1 = 0;
y5_solid1 = Y_out_core(1,1);
y6_solid1 = Y_out_core(2,1); 
y_solid1 = [y1_solid1; y2_solid1; y3_solid1; y4_solid1; y5_solid1; y6_solid1];

y1_solid2 = 1;
y2_solid2 = 0;
y3_solid2 = rho_l*g_out_core*y1_solid2;
y4_solid2 = 0;
y5_solid2 = 0;
y6_solid2 = 4*pi*G*rho_l*y1_solid2;
y_solid2 = [y1_solid2; y2_solid2; y3_solid2; y4_solid2; y5_solid2; y6_solid2];

y_solid3 = [0; 1; 0; 0; 0; 0];


Y_solid = [y_solid1, y_solid2, y_solid3];
Y_final = RK4_solid_incom(index_out+1, N, passo, r, n, G, rho_profile, ...
    mu_profile, Y_solid, matrix4,  r_in_core, r_out_core, ...
    r_mantle);


%SURFACE BOUNDARY CONDITIONS TO PRODUCE THE LINEAR COMBINATION OF y
%SOLUTIONS FOR THE CASE ---------------------------------------------------
g_N = get_g(a, r_in_core, r_out_core, r_mantle, rho_profile, G);

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

Cost = My \ CC;  %Cost = (C1, C2, C3) of the correct linear combination 
                 % at the surface


%LOVE NUMBERS------------------------------------------------------------
y1 = Cost(1)*Y_final(1,1) + Cost(2)*Y_final(1,2) + Cost(3)*Y_final(1,3);
y2 = Cost(1)*Y_final(2,1) + Cost(2)*Y_final(2,2) + Cost(3)*Y_final(2,3);
y5 = Cost(1)*Y_final(5,1) + Cost(2)*Y_final(5,2) + Cost(3)*Y_final(5,3);


switch Love_type
    case 'Shear'     %no potential, signs ok
        H = y1*g_N;
        L = y2*g_N;
        K = y5;
    otherwise        %potential with sign, signs to be changed
        H = -y1*g_N;
        L = -y2*g_N;
        K = y5 - 1;
end

% %Analytical LN for generic degree n, homogeneous incompressible case
% H_analitic = (2*n + 1)/(2*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));
% K_analitic = 3/(2*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));
% L_analitic = 3/(2*n*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));

end