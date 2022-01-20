function mgk=set_mgk(fi, lam, lam0 )
    fi=deg2rad(fi);
    lam=deg2rad(lam);
    lam0=deg2rad(lam0);
    a=6378137; 
    e2=0.00669437999013; 
    
    b2= a^2 * (1-e2);
    ep2= (a^2-b2)/b2 ;%  - e prim ^2

    t=tan (fi);
    n2= ep2 * (cos(fi))^2;
    deltalamb=lam-lam0;
    mgk= 1+ deltalamb^2/2 * cos(fi)^2 *(1+n2)+ deltalamb^4/24 * cos(fi)^4 * (5- 4 * t^2) ;

end
