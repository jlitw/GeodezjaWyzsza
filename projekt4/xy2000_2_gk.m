function [Xgk, Ygk]=xy2000_2_gk(X2000, Y2000)
    m02000= 0.999923;

    Xgk=X2000/m02000 ;
    if Y2000>= 8000000
        nr_strefy =8;
    elseif Y2000> 7000000
        nr_strefy =7;
    elseif Y2000>= 6000000
        nr_strefy =6;
    elseif Y2000>= 5000000
        nr_strefy =5;
    
    end

    Ygk=(Y2000- nr_strefy *1000000-500000)/m02000;
    [f, l ]= gk2fil(Xgk, Ygk, nr_strefy*3);
    [Xgk, Ygk]=to_gk(f, l, 19);

end