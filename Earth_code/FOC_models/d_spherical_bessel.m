function dj_n = d_spherical_bessel (x, n, j_n)
     j_nminus = spherical_bessel(x,n-1);

     dj_n = j_nminus - ((n+1)/x)*j_n;

end