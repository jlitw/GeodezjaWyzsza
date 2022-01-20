function [X_gk,Y_gk, mgk]=to_gk(fi, lam, lam0)%seems to be good
    fi = deg2rad(fi);
    lam= deg2rad(lam);

    lam0= deg2rad(lam0);

    a=6378137; 
    e2=0.00669437999013; 
    
    %1
    b2= a^2 * (1-e2);
    ep2= (a^2-b2)/b2 ;%  - e prim ^2

    %2
    deltalamb = lam-lam0;
    t=tan (fi);
    n2= ep2 * (cos(fi))^2;

    N=a/sqrt(1-(e2*(sin(fi))^2));
    
    %3 
    A0=1- e2/4-3 *e2^2/64-5 *e2^3/256 ;
    A2=3/8 * (e2+ e2^2/4+ 15* e2^3 /128) ;
    A4=15/256 * (e2^2 + 3* e2^3 /4);
    A6= 35* e2^3/3072;
    sigma= a* (A0 * fi-A2 *sin(2*fi) + A4* sin(4*fi) - A6*sin(6*fi) );
    
    %4
    X_gk= sigma + deltalamb^2 / 2 *N * sin(fi) * cos(fi) * (1+deltalamb^2 /12 * (cos(fi))^2 * ( 5-t^2 + 9 * n2 + 4 * n2^2)+ deltalamb^4/ 360 * (cos(fi))^4 * (61-58*t^2 + t^4+270* n2-330 *n2* t^2 ) );
    
    Y_gk= deltalamb *N* cos(fi) * (1+deltalamb^2 /6 * (cos(fi))^2 * ( 1- t^2 + n2)+ deltalamb^4/ 120 * (cos(fi))^4 * (5-18*t^2 + t^4 + 14*n2 -58 *n2* t^2 ) );

    mgk= 1+ deltalamb^2/2 * cos(fi)^2 *(1+n2)+ deltalamb^4/24 * cos(fi)^4 * (5- 4 * t^2) ;

end