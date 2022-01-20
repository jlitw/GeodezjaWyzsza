function [X_gk, Y_gk]=xy92_2_gk( X92, Y92)
    m092= 0.9993;
    x0=5300000.0;
    y0=500000.0;
    
    X_gk=(X92+x0)/m092 ;
    Y_gk=(Y92-y0)/m092 ;
end