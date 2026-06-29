function ddj_n = dd_spherical_bessel (x, n, j_n)
     j_nminus = spherical_bessel(x,n-1);

     ddj_n = (((n+1)*(n+2) - x^2)/(x^2))*j_n - (2/x)*j_nminus;

end