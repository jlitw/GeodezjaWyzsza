function [fiB, lamB, Aab]= Kivioj(fiA, lamA, s, A)
    a=6378137; 
    e2=0.00669437999013;
    fiA = deg2rad(fiA);
    lamA= deg2rad(lamA);
    A= deg2rad(A);

    n= round(s/1000);
    ds=s/n;
    
    for i=1:n
        N=a/sqrt(1-(e2*(sin(fiA))^2));
        M=(a *(1-e2)) / sqrt(( 1-(e2*(sin(fiA))^2))^3 );
        
        deltfi = ds * cos(A) / M;
        fi_m = fiA +0.5* deltfi ;
    
        A_m = A + 0.5*ds/N *sin(A) * tan(fi_m);
    
        N_m=a/sqrt(1-(e2*(sin(fi_m))^2));
        M_m=(a *(1-e2)) / sqrt(( 1-(e2*(sin(fi_m))^2))^3 );
    
        deltfi_m = ds *cos(A_m)/ M_m;
        fiA=fiA + deltfi_m;
    
        deltL_m = ds * sin(A_m) / (N_m * cos(fi_m));
        lamA = lamA +deltL_m;
    
        deltA_m = ds/N_m * sin(A_m) * tan(fi_m);
        A = A + deltA_m;
    
    end
    
    fiB= fiA;
    lamB = lamA;
    Aab = A;
    fiB= rad2deg(fiB);
    lamB =rad2deg(lamB);
    Aab =rad2deg (Aab);
end
