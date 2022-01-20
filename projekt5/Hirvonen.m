function [fi,lam, h] =Hirvonen (x, y, z )
    format longg
    %dane ->Krasowski
     e2 = 0.0066934215520398155;
     a = 6378245;
    % e2 = 0.00669342162296085;
    % a = 6378245;

    r = sqrt(x^2+y^2);
    
    eps = deg2rad(0.00005/3600);
    delt=1;
    n=0;
    while abs((delt))>eps 
        n=n+1;
        if n<=1
            fi1 = atan((z/r) * ((1 - e2)^-1));
            N1 = a/sqrt(1 - e2 * (sin(fi1^2)));
            h1=(r/cos(fi1))-N1;
        end
        if n>1
            fi2 = atan((z/r) * (1 - e2 * (N1/(N1 + h1)))^-1);
            N2 = a/sqrt(1 - e2 * (sin(fi2)^2));
            h2=(r/cos(fi2))-N2;
            delt=fi2-fi1;
            
            fi1=fi2;
            N1=N2;
            h1=h2;
        end
    end
    fi=fi2;
    lam=atan(y/x);
    N=a/ sqrt( 1-e2 * sin(fi)^2);
    h= r/cos(fi)-N ;

    lam= rad2deg(lam);
    fi=rad2deg(fi);
end