%-------------------------------------------------------------------------%
% --- Evaluate LNs of harmonic degree n and 'Love_type' type -------------%
% --- for compressible planet, in a range of the specified ---------------%
% --- elastic parameter --------------------------------------------------%
%-------------------------------------------------------------------------%

function[H, K, L] = main_compressible(n, Love_type, mu_mantle)

%PARAMETERS SETTINGS ------------------------------------------------------
G = 6.67 * 10^(-11);  %in [N m^2 / kg^2]

[a, r_core, r_mantle, rho_profile, mu_profile, lambda_profile] = ...
    set_parameters();

[passo, r0, r, N] = set_integration(a, r_core);

%select the elastic parameter to test in a range
mu_profile(2) = mu_mantle;


%VARIABLES SCALING --------------------------------------------------------
[a, G, r_core, r_mantle, rho_profile, mu_profile, lambda_profile, ...
    passo, r0, r] = normalize(a, G, r_core, r_mantle, rho_profile, ...
    mu_profile, lambda_profile, passo, r0, r);


%STARTING SOLUTIONS IN r0 -------------------------------------------------
rho_0 = get_parameter(r0, r_core, r_mantle, rho_profile);
mu_0 = get_parameter(r0, r_core, r_mantle, mu_profile);
lambda_0 = get_parameter(r0, r_core, r_mantle, lambda_profile);
g_0 = get_g(r0, r_core, r_mantle, rho_profile, G);
A_0 = lambda_0 + 2*mu_0;

epsilon = g_0/r0; %assumed to be constant in r0

v_S = sqrt(mu_0/rho_0);
v_P = sqrt(A_0/rho_0);
alpha2 = (4*pi*G*rho_0 + epsilon) / v_P^2;
gamma2 = (4*n*(n+1)*(epsilon^2)) / (v_P^2 * v_S^2);

k2 = 0.5*(alpha2 + sqrt(alpha2^2 + gamma2));
q2 = 0.5*(alpha2 - sqrt(alpha2^2 + gamma2));

%Starting solution 1
k = sqrt(k2);
kr = k*r0;
C_k = - epsilon / (v_S^2 * k^2);

j_n = spherical_bessel(kr, n);
dj_n = d_spherical_bessel(kr, n, j_n);
ddj_n = dd_spherical_bessel(kr, n, j_n);

y1_s0_1 = -(1/k) * (C_k*n*(n+1)*j_n/kr + dj_n);
y2_s0_1 = -(1/k) * ((1+C_k)*j_n/kr + C_k*dj_n); 
y3_s0_1 = lambda_0*j_n + 2*mu_0*C_k*n*(n+1)*(j_n/kr - dj_n)/kr - ...
    2*mu_0*ddj_n;
y4_s0_1 = -mu_0 * (C_k*j_n - 2*(1+C_k)*(j_n/kr - dj_n)/kr + ...
    2*C_k*ddj_n);
y5_s0_1 = 4*pi*G*rho_0*j_n / k^2;
y6_s0_1 = (1 - n*C_k)*4*pi*G*rho_0*(n+1)*j_n / (kr*k);

%Starting solution 2
if q2 < 0
    q = sqrt(-q2);
    qr = q * r0;
    i_n = spherical_bessel_modified(qr, n);
    di_n = d_spherical_bessel_modified(qr, n, i_n);
    ddi_n = dd_spherical_bessel_modified(qr, n, i_n);
else
    q = sqrt(q2);
    qr = q * r0;
    i_n = spherical_bessel(qr, n);
    di_n = d_spherical_bessel(qr, n, i_n);
    ddi_n = dd_spherical_bessel(qr, n, i_n);
end

C_q = - epsilon / (v_S^2 * q^2); 

y1_s0_2 = -(1/q) * (C_q*n*(n+1)*i_n/qr + di_n);
y2_s0_2 = -(1/q) * ((1+C_q)*i_n/qr + C_q*di_n); 
y3_s0_2 = lambda_0*i_n + 2*mu_0*C_q*n*(n+1)*(i_n/qr - di_n)/qr - ...
    2*mu_0*ddi_n;
y4_s0_2 = -mu_0 * (C_q*i_n - 2*(1+C_q)*(i_n/qr - di_n)/qr + ...
    2*C_q*ddi_n);
y5_s0_2 = 4*pi*G*rho_0*i_n / q^2;
y6_s0_2 = (1 - n*C_q)*4*pi*G*rho_0*(n+1)*i_n / (qr*q);


%Starting solution 3
y1_s0_3 = n * ((r0)^(n-1));
y2_s0_3 = r0^(n-1);
y3_s0_3 = 2*mu_0*n*(n-1)*((r0)^(n-2));
y4_s0_3 = 2*mu_0*(n-1)*((r0)^(n-2));
y5_s0_3 = -n*epsilon*((r0)^n);
y6_s0_3 = (4*pi*G*rho_0 - (2*n + 1)*epsilon)*n*((r0)^(n-1));

%column arrays of starting sets
Y_s0_1 = [y1_s0_1; y2_s0_1; y3_s0_1; y4_s0_1; y5_s0_1; y6_s0_1];
Y_s0_2 = [y1_s0_2; y2_s0_2; y3_s0_2; y4_s0_2; y5_s0_2; y6_s0_2];
Y_s0_3 = [y1_s0_3; y2_s0_3; y3_s0_3; y4_s0_3; y5_s0_3; y6_s0_3];

%matrix of starting sets
Y_s0 = [Y_s0_1, Y_s0_2, Y_s0_3];  


%NUMERICAL INTEGRATION VIA RUNGE KUTTA 4-----------------------------------
matrix4 = [0, 0, 0, 0, 0; 1/2, 1/2, 0, 0, 0; 1/2, 0, 1/2, 0, 0; ...
    1, 0, 0, 1, 0; 0, 1/6, 1/3, 1/3, 1/6];

Y_final = RK4_solid_com(2, N, passo, r, n, G, rho_profile, ...
    mu_profile, lambda_profile, Y_s0, matrix4, r_core, r_mantle);


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

end