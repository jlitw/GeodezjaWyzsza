function [P] = area (fi1, l1, fi2, l2)
    fi1 = deg2rad(fi1);
    l1= deg2rad(l1);

    fi2 = deg2rad(fi2);
    l2 = deg2rad(l2);

    a=6378137; 
    e2=0.00669437999013;
    b=a* sqrt(1-e2);
    e = sqrt(e2);

    wf1 = sin(fi1) /  (1 - e2 * sin(fi1)^2 ) + 1/ (2*e) * log( (1+e* sin(fi1))/(1-e* sin(fi1)));
    wf2 = sin(fi2) /  (1 - e2 * sin(fi2)^2 ) + 1/ (2*e) * log( (1+e* sin(fi2))/(1-e* sin(fi2)));

    P = b^2 * (l1-l2) /2 * (wf2 -wf1);

end