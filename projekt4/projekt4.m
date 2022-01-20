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

x=zeros(1,6);
y=zeros(1,6);
x92=zeros(1,6);
y92=zeros(1,6);
x2000=zeros(1,6);
y2000=zeros(1,6);

Xw=zeros(1,6);
Yw=zeros(1,6);
Zw=zeros(1,6);

fispr=zeros(1,6);
lamspr=zeros(1,6);


xsprn=zeros(1,6);
ysprn=zeros(1,6);
xn2=zeros(1,6);
yn2=zeros(1,6);

s=zeros(1,6);
mgk=zeros(1,6);
m92=zeros(1,6);
m2000=zeros(1,6);
kgk=zeros(1,6);
k92=zeros(1,6);
k2000=zeros(1,6);
kgk_2=zeros(1,6);
k92_2=zeros(1,6);
k2000_2=zeros(1,6);

n=1;
while n<7
    [s(n)]=strefa(l(n));
    [x(n), y(n), mgk(n)]= to_gk(f(n),l(n),19);
    [fispr(n), lamspr(n)]=gk2fil(x(n),y(n), 19);

    [x92(n),y92(n), m92(n)]= to_1992(x(n),y(n));
    [x2000(n),y2000(n), m2000(n)]= to_2000(x(n),y(n), s(n));


    [xsprn(n) ,ysprn(n)]=xy92_2_gk(x92(n), y92(n));
    [xn2(n), yn2(n)]= xy2000_2_gk(x2000(n), y2000(n) );
    
    kgk(n)= (mgk(n)-1)*1000;
    k92(n)=(m92(n)-1)*1000;
    k2000(n)=(m2000(n)-1)*1000;

    fprintf('\n %d', n );
    fprintf('\ndane phi: %f', f(n) );
    fprintf('\ndane lambda: %f', l(n));
    fprintf('\nGaus-K x: %f', x(n) );
    fprintf('\nGaus-K y: %f', y(n) );
  
    fprintf('\nSPR phi: %f', fispr(n) );
    fprintf('\nspr lambda: %f', lamspr(n));
    
    fprintf('\n1992 x: %f', x92(n) );
    fprintf('\n1992 y: %f', y92(n) );
    fprintf('\n2000 x: %f', x2000(n) );
    fprintf('\n2000 y: %f', y2000(n) );
     
    
    fprintf('\n1992-> gk x: %f', xsprn(n) );
    fprintf('\n1992 ->gky: %f', ysprn(n) );
    
    fprintf('\n2000-> gk x: %f', xn2(n) );
    fprintf('\n2000 ->gky: %f', yn2(n) );
    fprintf('\n\n ' );
    
    fprintf('\nMgk: %f', mgk(n) );
    fprintf('\nM1992: %f', m92(n));
    fprintf('\nM2000: %f', m2000(n));

    fprintf('\nKgk: %f', kgk(n) );
    fprintf('\nK92: %f', k92(n) );
    fprintf('\nK2000: %f', k2000(n) );
    fprintf('\n\n ' );
    fprintf('\nMgk^2: %f', mgk(n)^2 );
    fprintf('\nM1992^2: %f' , m92(n)^2);
    fprintf('\nM2000^2: %f' , m2000(n)^2);

    kgk_2(n)=(mgk(n)^2-1)*10000;
    k92_2(n)=(m92(n)^2-1)*10000;
    k2000_2(n)=(m2000(n)^2-1)*10000;

    fprintf('\nk^2: %f', kgk_2(n) );
    fprintf('\nk1992^2: %f' , k92_2(n));
    fprintf('\nk2000^2: %f' , k2000_2(n));

    fprintf('\n\n ' );

    n=n+1;
end

t=[x92(1),x92(3),x92(4),x92(2), x92(1)];
t2=[y92(1),y92(3),y92(4),y92(2),y92(1)];
plot(t,t2, '-');
axis equal
p92=polyarea(t, t2)/1000000;

t=[x(1),x(3),x(4),x(2), x(1)];
t2=[y(1),y(3),y(4),y(2),y(1)];
plot(t,t2, '-');
axis equal
pgk=polyarea(t, t2)/1000000;

t=[x2000(1),x2000(3),x2000(4),x2000(2), x2000(1)];
t2=[y2000(1),y2000(3),y2000(4),y2000(2),y2000(1)];
plot(t,t2, '-');
axis equal
p2000=polyarea(t, t2)/1000000;


pelipsoidalne= area(f(1),l(1), f(4),l(4))/1000000;

fprintf('\npole elipsoidalne: %f', pelipsoidalne );
fprintf(' km^2\npole gk %f', pgk);
fprintf(' km^2\npole 2000 %f', p2000 );
fprintf(' km^2\npole 1992: %f', p92 );
fprintf(' km^2 \n\n ' );

function [str]= strefa(lam)
 if(lam>=22.5)
     str=8;
 elseif(lam<22.5 && lam>=19.5)
     str=7;
 elseif(lam<19.5 && lam>=16.5)
     str=6;
 elseif(lam<16.5)
     str=5;
 end
end