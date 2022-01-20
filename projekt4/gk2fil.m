function [fi, lam] =gk2fil(X_gk, Y_gk , lam0)
    lam0=deg2rad(lam0);
    a=6378137; 
    e2=0.00669437999013; 
    fip=1000;
    
    A0=1-(e2/4)-((3*e2^2)/64)-((5*e2^3)/256);
    A2=3/8*(e2+(e2^2/4)+(15*e2^3/128));
    A4=15/256*(e2^2+(3*e2^3/4));
    A6=(35*e2^3)/3072;
    
    fi=X_gk/(a*A0);
    while abs(fi-fip)>0.000001*pi/(180*3600) 
        sigma=a*(A0*fi-A2*sin(2*fi)+A4*sin(4*fi)-A6*sin(6*fi));
        fip=fi;
        fi=fip+((X_gk-sigma)/(a*A0));
    end
    N=a/sqrt(1-(e2*(sin(fi))^2));
    M=(a *(1-e2)) / sqrt(( 1-(e2*(sin(fi))^2))^3 );
    t=tan(fi);
    b2=a^2*(1-e2);
    ep2=(a^2-b2)/b2;
    n2=ep2*cos(fi)^2;
    
    fik=fi-(Y_gk^2*t)/(2*M*N)*(1-(Y_gk^2)/(12*N^2)*(5+3*t^2+n2-9*n2*t^2-4*n2^2)+(Y_gk^4)/(360*N^4)*(61+90*t^2+45*t^4));
    lam=lam0+Y_gk/(N*cos(fi))*(1-(Y_gk^2/(6*N^2))*(1+2*t^2+n2)+(Y_gk^4)/(120*N^4)*(5+28*t^2+24*t^4+6*n2+8*n2*t^2));
    
    fi=rad2deg(fik) ;
    lam=rad2deg(lam) ;
end