function j_n = spherical_bessel (x, n)
     j_n = sqrt(pi/(2*x)) * besselj(n+0.5, x);
end