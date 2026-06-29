function i_n = spherical_bessel_modified(x, n)
    i_n = sqrt(pi/(2*x)) * besseli(n+0.5, x);
end