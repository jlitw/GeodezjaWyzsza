function [ sAB, Aab, Aba] = Vincent(fiA, lamA, fiB, lamB)
    a=6378137; 
    e2=0.00669437999013; 

    fiA = deg2rad(fiA);
    lamA= deg2rad(lamA);
    fiB = deg2rad(fiB);
    lamB = deg2rad(lamB);

    % 1
    b=a* sqrt(1-e2);
    f= 1- (b/a);
    
    %2
    deltlamb = lamB-lamA;
    
    %3
    Ua= atan((1-f)*tan(fiA));
    Ub= atan((1-f)*tan(fiB));
    
    %4
    L=deltlamb;
    Lb = 100;
    eps = 0.000001/3600;
    % 4 
    while abs(L-Lb)>=eps
        % 5
        sinsigma=(sqrt((cos(Ub)*sin(L))^2+(cos(Ua)*sin(Ub)-sin(Ua)*cos(Ub)*cos(L))^2));
        
        % 6
        cossigma= sin(Ua)*sin(Ub)+cos(Ua)*cos(Ub)*cos(L);
        
        % 7
        sigma= atan(sinsigma/cossigma);

        % 8
        sinalfa= cos(Ua)*cos(Ub)*sin(L)/ sinsigma ;
        
        % 9
        cos2alfa = 1 - sinalfa^2;
        
        % 10
        cos2sigmaM = cossigma - (2*sin(Ua)*sin(Ub) / (cos2alfa) );
        
        % 11
        C = f/16 *cos2alfa *  (4+ f *(4- 3* cos2alfa )) ;
        
        % 12
        Lb=L;
        L= deltlamb + (1-C) *f * sinalfa * ( sigma + C * sinsigma * ( cos2sigmaM + C * cossigma *(-1+2* (cos2sigmaM)^2 ) )   ) ;
    
    end

    % 13
    u2 = ((a^2 - b^2)* (cos2alfa)) / (b^2); 
    
    % 14
    A = 1 + u2/16384 * (4096 + u2 * (-768 + u2 * (320-175*u2) ) );

    % 15
    B = u2 / 1024 * ( 256 + u2 * ( -128 + u2 * ( 74- 47 * u2)));
    % 16
    deltasigma= B * sinsigma *(cos2sigmaM + 0.25* B * ( cossigma *(-1 + 2*(cos2sigmaM)^2 ) - (1/6 * B) * cos2sigmaM * (-3+4*(sinsigma)^2 )*(-3 + 4* (cos2sigmaM)^2) ));

    % 17
    sAB = b*A*( sigma - deltasigma); 

    % 18 19
    Aab = atan2 ( (cos(Ub) * sin(L) ), (cos(Ua) * sin(Ub) - sin(Ua) * cos(Ub) * cos(L) )  );
    Aba = atan2 ( (cos(Ua) * sin(L) ), (-1*sin(Ua) * cos(Ub) + cos(Ua) * sin(Ub) * cos(L) )) + pi;
    Aab =rad2deg(Aab);
    Aba= rad2deg(Aba);

    g=(cos(Ua) * sin(L) );
    d=(-1*sin(Ua) * cos(Ub) + cos(Ua) * sin(Ub) * cos(L) );
    if (g>0 && d>0 )
        Aba = atand (g / d) + 180;
    elseif(g>0 && d<0 )
        Aba =atand (g / d) + 360;
    elseif(g<0 && d<0 )
        Aba =atand (g / d) + 180;
    elseif(g<0 && d>0 )
        Aba =atand (g / d) + 540;
    end

    g2=(cos(Ub) * sin(L) );
    d2= (cos(Ua) * sin(Ub) - sin(Ua) * cos(Ub) * cos(L) )  ;
    if (g2>0 && d2>0 )
        Aab = atand (g2 / d2) ;
    elseif(g2>0 && d2<0 )
        Aab =atand (g2 / d2)+180 ;
    elseif(g2<0 && d2<0 )
        Aab =atand (g2 / d2) ;
    elseif(g2<0 && d2>0 )
        Aab =atand (g2 / d2) +360;
    end


end