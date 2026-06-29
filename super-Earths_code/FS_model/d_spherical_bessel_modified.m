%-------------------------------------------------------------------------%
% --- Produce the first derivative of a modified spherical Bessel --------%
% --- function i_n in x variable, of degree n ----------------------------%
%-------------------------------------------------------------------------%

function di_n = d_spherical_bessel_modified(x, n, i_n)
    in_minus = spherical_bessel_modified(x, n-1);

    di_n = in_minus - ((n+1)/x)*i_n;
end