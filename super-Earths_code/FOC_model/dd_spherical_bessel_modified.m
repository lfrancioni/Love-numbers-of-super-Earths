%-------------------------------------------------------------------------%
% --- Produce the second derivative of a modified spherical Bessel -------%
% --- function j_n in x variable, of degree n ----------------------------%
%-------------------------------------------------------------------------%

function ddi_n = dd_spherical_bessel_modified(x, n, i_n)
    in_minus = spherical_bessel_modified(x, n-1);

    ddi_n = (((n+1)*(n+2) + x^2)/(x^2))*i_n - (2/x)*in_minus;
end