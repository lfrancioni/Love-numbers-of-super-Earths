function ddj_n = dd_spherical_bessel (x, n, j_n)
     j_nminus = spherical_bessel(x,n-1);

     ddj_n = (((n*(n+1) + 2*(n+1)) / x^2) - 1)*j_n - (2/x)*j_nminus;

end