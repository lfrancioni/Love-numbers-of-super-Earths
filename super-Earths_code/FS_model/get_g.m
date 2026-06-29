%-------------------------------------------------------------------------%
% --- Calculate the gravity acceleration at radial distance r inside -----%
% --- the differentiated planet ------------------------------------------%
%-------------------------------------------------------------------------%

function g = get_g(r, r_core, r_mantle, rho_profile, G)
    if r < r_core
        M = rho_profile(1)*(4/3)*pi*r^3;
    elseif r < r_mantle
        M = rho_profile(1)*(4/3)*pi*(r_core)^3 + rho_profile(2)*(4/3)*pi* ...
            (r^3 - r_core^3);
    else
        M = rho_profile(1)*(4/3)*pi*(r_core)^3 + rho_profile(2)*(4/3)*pi*...
            (r_mantle^3 - r_core^3) + rho_profile(3)*(4/3)*pi*(r^3 - ...
            r_mantle^3);
    end

    g = (G*M)/r^2;
    
end
