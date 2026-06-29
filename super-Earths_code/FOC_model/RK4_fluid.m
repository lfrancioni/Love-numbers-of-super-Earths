%-------------------------------------------------------------------------%
% --- Solve the 4th order Runge-Kutta method for fluid -------------------%
% --- Y solution in each step of the interval [index_i, index_f] ---------%
%-------------------------------------------------------------------------%

function [dY] = RK4_fluid(index_i, index_f, step, r, n, G, ...
    rho, rho_profile, Y_s, butcher, r_in_core, r_out_core, r_mantle)

Y = Y_s;

for i = index_i:index_f 
        r1 = r(i-1);
        y5_step1 = Y(1, 1);   
        y7_step1 = Y(2, 1);   

        y_step1 = [y5_step1, y7_step1];
        f_step1 = system_fluid(r1, y_step1, n, G, rho, rho_profile, ...
            r_in_core, r_out_core, r_mantle);
        %row array: (f5_step1, f7_step1)


        r2 = r(i-1) + step * butcher(2,1);
        y5_step2 = y5_step1 + step * butcher(2,2) * f_step1(1);
        y7_step2 = y7_step1 + step * butcher(2,2) * f_step1(2);
   
        y_step2 = [y5_step2, y7_step2];
        f_step2 = system_fluid(r2, y_step2, n, G, rho, rho_profile, ...
            r_in_core, r_out_core, r_mantle);


        r3 = r(i-1) + step * butcher(3,1);
        y5_step3 = y5_step1 + step * (butcher(3,2) * f_step1(1) + ...
            butcher(3,3) * f_step2(1));
        y7_step3 = y7_step1 + step * (butcher(3,2) * f_step1(2) + ...
            butcher(3,3) * f_step2(2));

        y_step3 = [y5_step3, y7_step3];
        f_step3 = system_fluid(r3, y_step3, n, G, rho, rho_profile, ...
            r_in_core, r_out_core, r_mantle);


        r4 = r(i-1) + step * butcher(4,1);
        y5_step4 = y5_step1 + step * (butcher(4,2) * f_step1(1) + ...
            butcher(4,3) * f_step2(1) + butcher(4,4) * f_step3(1));
        y7_step4 = y7_step1 + step * (butcher(4,2) * f_step1(2) + ...
            butcher(4,3) * f_step2(2) + butcher(4,4) * f_step3(2));

        y_step4 = [y5_step4, y7_step4];
        f_step4 = system_fluid(r4, y_step4, n, G, rho, rho_profile, ...
            r_in_core, r_out_core, r_mantle);

        y_RK = zeros(1, 2);  %solution array of num. integration
        for m = 1:2
            y_RK(m) = Y(m, 1) + step * (butcher(5,2) * f_step1(m) + butcher(5,3) * f_step2(m) + ...
                butcher(5,4) * f_step3(m) + butcher(5,5) * f_step4(m));
            Y(m, 1) = y_RK(m);
        end

end

dY = Y;

end
