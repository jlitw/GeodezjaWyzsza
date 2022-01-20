function  [X92, Y92, mu]= to_1992(X_gk, Y_gk)
    [fi, lam]=gk2fil(X_gk, Y_gk,19);
    m092= 0.9993;
    x0=5300000.0;
    y0=500000.0;
    
    X92=m092* X_gk - x0;
    Y92=m092*Y_gk + y0;

    mgk=set_mgk(fi, lam, 19);
    mu=mgk*m092;

end
