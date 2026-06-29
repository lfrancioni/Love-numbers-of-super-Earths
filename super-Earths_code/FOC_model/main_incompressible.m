%-------------------------------------------------------------------------%
% --- Evaluate LNs of harmonic degree n and 'Love_type' type -------------%
% --- for incompressible planet, in a range of the specified -------------%
% --- elastic parameter --------------------------------------------------%
%-------------------------------------------------------------------------%

function [H, K, L] = main_incompressible(n, Love_type, mu_mantle)

%PARAMETERS SETTINGS-------------------------------------------------------
G = 6.67 * 10^(-11);  %in [N m^2 / kg^2]

[a, r_in_core, r_out_core, r_mantle, rho_profile, mu_profile, lambda_profile] = ...
    set_parameters();
%lambda_profile will be unused in this incompressible case

[step, r0, r, N] = set_integration(a);

%select the elastic parameter to test in a range
mu_profile(3) = mu_mantle;

%indexes for num. integration
index_in = find(r == r_in_core);   
index_out = find(r == r_out_core);


%VARIABLES SCALING --------------------------------------------------------
%(inspired by: prepare_planet_model.py - LoadDef by H.Martens)
[a, G, r_in_core, r_out_core, r_mantle, rho_profile, mu_profile, ...
    lambda_profile, step, r0, r] = normalize(a, G, r_in_core, ...
    r_out_core, r_mantle, rho_profile, mu_profile, lambda_profile, ...
    step, r0, r);


%STARTING SOLUTIONS IN r0--------------------------------------------------
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

%Column arrays of starting sets
Y_s0_1 = [y1_s0_1; y2_s0_1; y3_s0_1; y4_s0_1; y5_s0_1; y6_s0_1];
Y_s0_2 = [y1_s0_2; y2_s0_2; y3_s0_2; y4_s0_2; y5_s0_2; y6_s0_2];
Y_s0_3 = [y1_s0_3; y2_s0_3; y3_s0_3; y4_s0_3; y5_s0_3; y6_s0_3];

%matrix of starting sets
Y_s0 = [Y_s0_1, Y_s0_2, Y_s0_3];


%NUMERICAL INTEGRATION VIA RUNGE KUTTA 4 ----------------------------------
matrix4 = [0, 0, 0, 0, 0; 1/2, 1/2, 0, 0, 0; 1/2, 0, 1/2, 0, 0; 1, 0, 0, 1, 0; ...
    0, 1/6, 1/3, 1/3, 1/6];

%Inner core integration
Y_in_core = RK4_solid_incom(2, index_in, step, r, n, G, rho_profile, ...
    mu_profile, Y_s0, matrix4, r_in_core, r_out_core, r_mantle);

%Outer core integration
%B.C. entering a fluid layer
rho_l = rho_profile(2);
g_in_core = get_g(r_in_core, r_in_core, r_out_core, r_mantle, rho_profile, G);

% In incompressible case, the b.cs. employed for the compressible case
% are singular due to: Y_in_core(4,3)=0 , so...

%Null space of shear stress y_4 combination
A_shear = [Y_in_core(4,1), Y_in_core(4,2), Y_in_core(4,3)];
V = null(A_shear); % matrix 3x2

%Combining the 3 null combinations in 2 shear-free Z solutions
Z = Y_in_core * V; 

%W constant evaluation from the Z solutions
Wn = Z(1,1) + Z(5,1)/g_in_core - Z(3,1)/(rho_l * g_in_core);
Wd = Z(1,2) + Z(5,2)/g_in_core - Z(3,2)/(rho_l * g_in_core);

W = Wn / Wd;

%Variables at the solid-fluid interface
y3_0 = Z(3,1) - W * Z(3,2);
y5_0 = Z(5,1) - W * Z(5,2);
y6_0 = Z(6,1) - W * Z(6,2);

y7_0 = y6_0 - 4*pi*G*y3_0/g_in_core;

Y_fluid = [y5_0; y7_0];
Y_out_core = RK4_fluid(index_in+1, index_out, step, r, n, G, ...
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
Y_final = RK4_solid_incom(index_out+1, N, step, r, n, G, rho_profile, ...
    mu_profile, Y_solid, matrix4,  r_in_core, r_out_core, ...
    r_mantle);


%SURFACE BOUNDARY CONDITIONS TO PRODUCE THE LINEAR COMBINATION OF y
%SOLUTIONS ----------------------------------------------------------------
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

Const = My \ CC;  %Const = (C1, C2, C3) of the correct linear combination 
                  %at the surface


%LOVE NUMBERS------------------------------------------------------------
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

% %Analytical LN for generic degree n, homogeneous incompressible case
% H_analitic = (2*n + 1)/(2*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));
% K_analitic = 3/(2*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));
% L_analitic = 3/(2*n*(n-1)*(1 + (2*n^2 + 4*n + 3)*mu_profile(1)/...
%     (n*rho_profile(1)*g(N)*a)));

end