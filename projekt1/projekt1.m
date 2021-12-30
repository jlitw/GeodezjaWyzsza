clear; 
a=6378137; 
e2=0.00669437999013; 

%wsp samolotu 
macierzDane=load('dane.txt');
phi=macierzDane(:,1);  
lambda=macierzDane(:,2); 
h=macierzDane(:,3); 

%wsp lotniska 
%dane   33.9361	-118.4166	84
phiL=33.9361;
lambdaL=-118.4166;
hL=84;									


%wsp.lotnisko 
N_L=a/sqrt(1-(e2*(sind(phiL))^2));
xL=(N_L+hL)*cosd(phiL)*cosd(lambdaL);
yL=(N_L+hL)*cosd(phiL)*sind(lambdaL);
zL=(N_L*(1-e2)+hL)*sind(phiL);

%wsp. samototu xyz
N=(a./(1-(e2*(sind(phi)).^2)).^(0.5) );
x=(N+h).*cosd(phi).*cosd(lambda);
y=(N+h).*cosd(phi).*sind(lambda);
z=(N*(1-e2)+h).*sind(phi);


%wsp neu
m_delt=[x-xL y-yL z-zL]';
m_tr=[- sind(phiL)*cosd(lambdaL) -sind(lambdaL) cosd(phiL)*cosd(lambdaL)
    - sind(phiL)*sind(lambdaL) cosd(lambdaL) cosd(phiL)*sind(lambdaL)
    cosd(phiL) 0    sin(phiL)]';

macierz_neu=m_tr*m_delt;

%wspolrzedne neu 
n=macierz_neu(1,:);
e=macierz_neu(2,:);
u=macierz_neu(3,:);

%s- odleglosc skosna
s=(n.^2+e.^2+u.^2).^(0.5);

%A-azymut
A=atand(e./n);
N=562;
for j=1:N
     if A(1,j)<0
        A(1,j)=A(1,j)+360;
     end
end



%z - odleglosc zenitalna obiektu 
z=acosd(u./(n.^2+e.^2+u.^2).^(0.5));


%wykresy 
nexttile
geoscatter(phi,lambda,5, 'ro')

nexttile
plot3(n,e,u);
grid;
title('Wykres krzywej w przestrzeni trójwymiarowej')
xlabel('n');
ylabel('e');
zlabel('u') ;
 


nexttile
plot(s, u);
grid;
title('Wykres s u')
xlabel('s');
ylabel('u');




