function [x, y,z] = fl_to_xyz(phi, lambda, h)
    format longg
    %grs80
    e2 = 0.00669437999013;
    a = 6378137;

    phi = deg2rad(phi);
    lambda = deg2rad(lambda);

    N=a/sqrt(1-e2*(sin(phi))^2);
    x=(N+h)*cos(phi)*cos(lambda);
    y=(N+h)*cos(phi)*sin(lambda);
    z=(N*(1-e2)+h)*sin(phi);

end