function [Xw, Yw, Zw]= tr_Bursy_Wolfa(X, Y, Z)
    format longg

    x0=-33.4297;
    y0=146.5746;
    z0=76.2865;
    m=1+0.8407728/10^6;
    ex=deg2rad(-0.35867/3600);
    ey=deg2rad(-0.05283/3600);
    ez=deg2rad(0.84354/3600);

    m_poczat=[X; 
        Y; 
        Z];
    m = [m ez -ey;
        -ez m  ex;
        ey -ex m];
    m0= [x0; 
        y0; 
        z0];

    wynik = m*m_poczat + m0;
    Xw=wynik(1);
    Yw=wynik(2);
    Zw=wynik(3);

end