%-------------------------------------------------------------------------%
% --- Solve the 4th order Runge-Kutta method for compressible ------------%
% --- Y solution in each step of the interval [index_i, index_f] ---------%
%-------------------------------------------------------------------------%

function [dY] = RK4_solid_com(index_i, index_f, step, r, n, G, rho_profile, ...
    mu_profile, lambda_profile, Y_s, butcher, r_core, r_mantle)

Y = Y_s;

for i = index_i:index_f
    r1 = r(i-1);
    %Use always the same set of elastic parameters in each Runge-Kutta
    %step ("freezing"), not adapting them step by step
    rho = get_parameter(r1, r_core, r_mantle, rho_profile);
    mu = get_parameter(r1, r_core, r_mantle, mu_profile);
    lambda = get_parameter(r1, r_core, r_mantle, lambda_profile);
    

    for l = 1:3
        y1_step1 = Y(1, l);     %y1 = u
        y2_step1 = Y(2, l);     %y2 = v
        y3_step1 = Y(3, l);     %y3 = T_r
        y4_step1 = Y(4, l);     %y4 = T_theta
        y5_step1 = Y(5, l);     %y5 = Phi
        y6_step1 = Y(6, l);     %y6 = Q

        y_step1 = [y1_step1, y2_step1, y3_step1, y4_step1, y5_step1, y6_step1];
        f_step1 = system_solid_com(r1, y_step1, n, G, rho, rho_profile, ...
            mu, lambda, r_core, r_mantle);
        %row array: (f1_step1, f2_step1, ... , f6_step1)


        r2 = r1 + step * butcher(2,1);
        y1_step2 = y1_step1 + step * butcher(2,2) * f_step1(1);
        y2_step2 = y2_step1 + step * butcher(2,2) * f_step1(2);
        y3_step2 = y3_step1 + step * butcher(2,2) * f_step1(3);
        y4_step2 = y4_step1 + step * butcher(2,2) * f_step1(4);
        y5_step2 = y5_step1 + step * butcher(2,2) * f_step1(5);
        y6_step2 = y6_step1 + step * butcher(2,2) * f_step1(6);

        y_step2 = [y1_step2, y2_step2, y3_step2, y4_step2, y5_step2, y6_step2];
        f_step2 = system_solid_com(r2, y_step2, n, G, rho, rho_profile, ...
            mu, lambda, r_core, r_mantle);


        r3 = r1 + step * butcher(3,1);
        y1_step3 = y1_step1 + step * (butcher(3,2) * f_step1(1) + ...
            butcher(3,3) * f_step2(1));
        y2_step3 = y2_step1 + step * (butcher(3,2) * f_step1(2) + ...
            butcher(3,3) * f_step2(2));
        y3_step3 = y3_step1 + step * (butcher(3,2) * f_step1(3) + ...
            butcher(3,3) * f_step2(3));
        y4_step3 = y4_step1 + step * (butcher(3,2) * f_step1(4) + ...
            butcher(3,3) * f_step2(4));
        y5_step3 = y5_step1 + step * (butcher(3,2) * f_step1(5) + ...
            butcher(3,3) * f_step2(5));
        y6_step3 = y6_step1 + step * (butcher(3,2) * f_step1(6) + ...
            butcher(3,3) * f_step2(6));

        y_step3 = [y1_step3, y2_step3, y3_step3, y4_step3, y5_step3, y6_step3];
        f_step3 = system_solid_com(r3, y_step3, n, G, rho, rho_profile, ...
            mu, lambda, r_core, r_mantle);


        r4 = r1 + step * butcher(4,1);
        y1_step4 = y1_step1 + step * (butcher(4,2) * f_step1(1) + ...
            butcher(4,3) * f_step2(1) + butcher(4,4) * f_step3(1));
        y2_step4 = y2_step1 + step * (butcher(4,2) * f_step1(2) + ...
            butcher(4,3) * f_step2(2) + butcher(4,4) * f_step3(2));
        y3_step4 = y3_step1 + step * (butcher(4,2) * f_step1(3) + ...
            butcher(4,3) * f_step2(3) + butcher(4,4) * f_step3(3));
        y4_step4 = y4_step1 + step * (butcher(4,2) * f_step1(4) + ...
            butcher(4,3) * f_step2(4) + butcher(4,4) * f_step3(4));
        y5_step4 = y5_step1 + step * (butcher(4,2) * f_step1(5) + ...
            butcher(4,3) * f_step2(5) + butcher(4,4) * f_step3(5));
        y6_step4 = y6_step1 + step * (butcher(4,2) * f_step1(6) + ...
            butcher(4,3) * f_step2(6) + butcher(4,4) * f_step3(6));

        y_step4 = [y1_step4, y2_step4, y3_step4, y4_step4, y5_step4, y6_step4];
        f_step4 = system_solid_com(r4, y_step4, n, G, rho, rho_profile, ...
            mu, lambda, r_core, r_mantle);

        
        y_RK = zeros(1, 6);  %solution array of num. integration
        for m = 1:6
            y_RK(m) = Y(m, l) + step * (butcher(5,2) * f_step1(m) + ...
                butcher(5,3) * f_step2(m) + butcher(5,4) * f_step3(m) + ...
                butcher(5,5) * f_step4(m));
            Y(m, l) = y_RK(m);
        end

    end

end

dY = Y;

end
