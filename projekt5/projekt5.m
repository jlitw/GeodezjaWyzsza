clear;
format longg

clear ;
f(1)=50.25;
l(1)=20.75;
f(2)=50;
l(2)=20.75;
f(3)=50.25;
l(3) =21.25;
f(4)=50;
l(4) =21.25;
f(5)=50.125000;
l(5) =21.000000;
f(6)=50.125269;
l(6)=21.000675;
h=100;

x=zeros(1,6);
y=zeros(1,6);
z=zeros(1,6);

Xw=zeros(1,6);
Yw=zeros(1,6);
Zw=zeros(1,6);

fi=zeros(1,6);
lam=zeros(1,6);
h2=zeros(1,6);

xS=zeros(1,6);
yS=zeros(1,6);
zS=zeros(1,6);
N=zeros(1,6);

e2 = 0.0066934215520398155;
a = 6378245;


n=1;
while n<7
    [x(n),y(n), z(n)]= fl_to_xyz(f(n), l(n), h);
    [Xw(n), Yw(n), Zw(n)]=tr_Bursy_Wolfa(x(n), y(n),z(n));
    [fi(n),lam(n),h2(n)]= Hirvonen(Xw(n), Yw(n),Zw(n));
    N(n)=a/sqrt(1-e2*(sind(fi(n)))^2);

    fprintf('\n %d', n );
    fprintf('\ndane phi: %f', f(n) );
    fprintf('\ndane lambda: %f', l(n));
    fprintf('\ndane h: %f', h);
    fprintf('\n GRS80 w xyz: ' );
    fprintf('\nx: %9.3f', x(n) );
    fprintf('\ny: %9.3f', y(n) );
    fprintf('\nz: %9.3f', z(n));

    fprintf('\n Krasowski: \nx: %9.3f', Xw(n) );
    fprintf('\ny: %9.3f', Yw(n) );
    fprintf('\nz: %9.3f', Zw(n) );
  
    fprintf('\nK phi: %f', fi(n) );
    fprintf('\nK lambda: %f', lam(n));
    fprintf('\nK h: %f', h2(n));

    xS(n)=(N(n)+h2(n))*cosd(fi(n))*cosd(lam(n));
    yS(n)=(N(n)+h2(n))*cosd(fi(n))*sind(lam(n));
    zS(n)=(N(n)*(1-e2)+h2(n))*sind(fi(n));
    fprintf('\n Spr Hirvonen: \nx: %9.3f', xS(n) );
    fprintf('\ny: %9.3f', yS(n) );
    fprintf('\nz: %9.3f', zS(n) );

    fprintf('\n\n ' );

    n=n+1;
end
