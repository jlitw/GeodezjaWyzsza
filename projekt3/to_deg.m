function [st, m, sec]= to_deg(kat)
    st= floor(kat);
    kat=kat-st;
    m= floor(kat*60);
    kat=kat-m/60;
    sec=(kat*3600);
end

