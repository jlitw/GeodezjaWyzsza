clear; 
v=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23];

%dane gwiazdy rektastenzja(alfa) i deklinacja (delta) Zubeneschamali
alf = 15.284408333;%w h   15h 17m 03.87s 
delt = -9.386444445; %w stopniach  -9 23 11,2

% miejsca na Ziemi - dane
phi1=52.229769;  
lambda1=21.000196;

%Astorga | S 23° 13' 57.00" | W 51° 39' 56.02"
phi2= -23.2325;
lambda2= -51.665561111;

%Kisangani | N 0° 30' 55.01" | E 25° 11' 27.57"
phi3=0.515280555;
lambda3=25.190991667;

x_1=zeros(1,24);
y_1=zeros(1,24);
z_1=zeros(1,24);

x_2=zeros(1,24);
y_2=zeros(1,24);
z_2=zeros(1,24);

x_3=zeros(1,24);
y_3=zeros(1,24);
z_3=zeros(1,24);

z1=zeros(1,24); %z - odległosc zenitalna 
z2=zeros(1,24);
z3=zeros(1,24);

A1=zeros(1,24);
A2=zeros(1,24);
A3=zeros(1,24);

for h = 0:23
    t1 = katgodz(2001,19,10,h,lambda1,alf);
    t2 = katgodz(2001,19,10,h,lambda2,alf);
    t3 = katgodz(2001,19,10,h,lambda3,alf);
   
    %współrzedne horyzontalne:
    z1(h+1) = acosd(sind(phi1)*sind(delt)+cosd(phi1)*cosd(delt)*cosd(t1));
    A1(h+1) = atan2d ((-cosd(delt)*sind(t1)),(cosd(phi1)*sind(delt)-sind(phi1)*cosd(delt)*cosd(t1)));
    A1(h+1)= azymut(A1(h+1));

    z2(h+1) = acosd(sind(phi2)*sind(delt)+cosd(phi2)*cosd(delt)*cosd(t2));
    A2(h+1) = atan2d((-cosd(delt)*sind(t2)),(cosd(phi2)*sind(delt)-sind(phi2)*cosd(delt)*cosd(t2)));
    A2(h+1)= azymut(A2(h+1));

    z3(h+1) = acosd(sind(phi3)*sind(delt)+cosd(phi3)*cosd(delt)*cosd(t3));
    A3(h+1) = atan2d((-cosd(delt)*sind(t3)),(cosd(phi3)*sind(delt)-sind(phi3)*cosd(delt)*cosd(t3)));
    A3(h+1)= azymut(A3(h+1));

    r=1;
    %wspólrzedne xyz
    x_1(h+1) = r*sind(z1(h+1))*cosd(A1(h+1));
    y_1(h+1) = r*sind(z1(h+1))*sind(A1(h+1));
    z_1(h+1) = r*cosd(z1(h+1));

    x_2(h+1) = r*sind(z2(h+1))*cosd(A2(h+1));
    y_2(h+1) = r*sind(z2(h+1))*sind(A2(h+1));
    z_2(h+1) = r*cosd(z2(h+1));

    x_3(h+1) = r*sind(z3(h+1))*cosd(A3(h+1));
    y_3(h+1) = r*sind(z3(h+1))*sind(A3(h+1));
    z_3(h+1) = r*cosd(z3(h+1));
end

teta=linspace(0, 180, 19);
phi=linspace(0, 360, 36);
[t, p]=meshgrid(teta,phi);%
R=ones(size (t));

x=R.*cosd(t).*cosd(p);
y=R.*cosd(t).*sind(p);
z=R.*sind(t);

%wykresy 1
figure;
mesh(x,y,z)
axis equal
hold on;
plot3(x_1,y_1,z_1,'rx');
title('Wykres 1 - Warszawa');
xlabel('x');
ylabel('y');
zlabel('z');

figure;
nexttile
plot(v,90-z1,'-');
title('Wykres wysokości od czasu 1 - Warszawa')
xlabel('czas [h]');
ylabel('wysokość [°]');
grid;
nexttile
plot(v,A1,'-');
title('Wykres azymutu od czasu 1 - Warszawa')
xlabel('czas [h]');
ylabel('azymut [°]');
grid;

%wykresy 2
figure;
mesh(x,y,z)
axis equal
hold on;
plot3(x_2,y_2,z_2,'rx');
title('Wykres 2 - Astorga')
xlabel('x');
ylabel('y');
zlabel('z');

figure;
nexttile
plot(v,90-z2,'-');
title('Wykres wysokości od czasu 2 - Astorga')
xlabel('czas [h]');
ylabel('wysokość [°]');
grid;
nexttile
plot(v,A2,'-');
title('Wykres azymutu od czasu 2 - Astorga')
xlabel('czas [h]');
ylabel('azymut [°]');
grid;

%wykresy 3
figure;
mesh(x,y,z)
axis equal
hold on;
plot3(x_3,y_3,z_3,'rx');
title('Wykres 3 - Kisangani')
xlabel('x');
ylabel('y');
zlabel('z');

figure;
nexttile
plot(v, 90-z3,'-');
title('Wykres wysokości od czasu 3 - Kisangani')
xlabel('czas [h]');
ylabel('wysokość [°]');
grid;

nexttile
plot(v, A3,'-');
title('Wykres azymutu od czasu 3 - Kisangani')
xlabel('czas [h]');
ylabel('azymut [°]');
grid;

function [t] = katgodz(y,d,m,h,lambda,alfa) 
    jd=juliandate(datetime(y, m, d));
    g = GMST(jd); %stopnie 
    UT1 = h*1.002737909350795; %godziny 
    %obliczenie czasu gwiazdowego(w stopniach) 
    S = UT1*15 + lambda + g; 
    %obliczenie kąta godzinnego(w stopniach) 
    t = S - alfa*15; 
end  

function g= GMST(JD)
    T=(JD - 2451545)/36525;
    g = 280.46061837 + 360.98564736629 * (JD - 2451545) + 0.000387933 * T.^2 - T.^3/38710000;
    g=mod(g, 360);
end

function [x] = azymut(x)
    if x<0
        x=x+360;
    end
end