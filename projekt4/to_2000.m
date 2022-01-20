function  [X2000, Y2000, mu]=to_2000(X_gk, Y_gk, nr_strefy )
    [fi, lam]=gk2fil(X_gk, Y_gk,19);
    [x, y]=to_gk(fi, lam,nr_strefy*3);

    m02000= 0.999923;
    X2000=m02000* x;
    Y2000=m02000*y +nr_strefy *1000000 + 500000;

    mgk=set_mgk(fi, lam, 21);
    mu=mgk*m02000;

end
