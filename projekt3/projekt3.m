clear; 
a=6378137; 
e2=0.00669437999013; 

fA=51.25;
lA=20.75;
fB=51;
lB=20.75;
fC=51.25;
lC=21.25;
fD=51;
lD=21.25;

% Wyznaczyć punkt średniej szerokości 
fpSS=(fA+fD)/2;
lpSS=(lA+lD)/2;

% Wyznaczyć punkt środkowy przy użyciu algorytmu Vincentego i Kivioji 
% Wyznaczyć różnicę odległości pomiędzy tymi punktami 
% Wyznaczyć azymuty w tych punktach. 
[sAD, Aad, Ada]= Vincent(fA, lA, fD, lD);
[fi_K, l_K, Ak]=Kivioj(fA, lA, sAD/2, Aad);
[odleglosc, Ask, Aks] = Vincent(fi_K, l_K, fpSS, lpSS);


% Obliczyć pole powierzchni tego czworokąta 
[p]= area(fA, lA, fD, lD);

[fpSS(1),fpSS(2),fpSS(3)]=to_deg(fpSS);
fprintf('\nPunkt średniej szerokości (s) phi: %d', fpSS(1));
fprintf('° %d', fpSS(2));
fprintf("' %2.1f",fpSS(3));
[lpSS(1),lpSS(2),lpSS(3)]=to_deg(lpSS);
fprintf('"\t lambda: %d', lpSS(1));
fprintf('° %d', lpSS(2));
fprintf("' %2.1f",lpSS(3));


[Aad(1),Aad(2),Aad(3)]=to_deg(Aad);
fprintf('"\nAzymut AD: %d', Aad(1));
fprintf('° %d', Aad(2));
fprintf("' %2.5f",Aad(3));
[Ada(1),Ada(2),Ada(3)]=to_deg(Ada);
fprintf('"\nAzymut DA: %d', Ada(1));
fprintf('° %d', Ada(2));
fprintf("' %2.5f",Ada(3));


[fi_K(1),fi_K(2),fi_K(3)]=to_deg(fi_K);
fprintf('"\nPunkt środkowy (k) phi : %d', fi_K(1));
fprintf('° %d', fi_K(2));
fprintf("' %2.1f",fi_K(3));
[l_K(1),l_K(2),l_K(3)]=to_deg(l_K);
fprintf('"\t lambda: %d', l_K(1));
fprintf('° %d', l_K(2));
fprintf("' %2.1f",l_K(3));

fprintf('"\nOdległość pomiędzy p. średniej szerokości, a p. środkowym: %10.3f', odleglosc);

[Ask(1),Ask(2),Ask(3)]=to_deg(Ask);
fprintf(' m \nAzymut sk: %d', Ask(1));
fprintf('° %d', Ask(2));
fprintf("' %2.5f",Ask(3));
[Aks(1),Aks(2),Aks(3)]=to_deg(Aks);
fprintf('"\nAzymut ks: %d', Aks(1));
fprintf('° %d', Aks(2));
fprintf("' %2.5f",Aks(3));

fprintf('"\nPole powierzchni: %f', p);
fprintf(' m^2\n');
