clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
fontSize = 22;

v = [0 0 1; 
     0 1 0; 
     1 0 0]; 
 
v_rot = [0 0 0];
v_proj = [0 0];
v_proj_KS = [0 0];
v_proj_temp = [0 0];
Varianten = [0 0 0];
Bild = [0 0];
 
Messdaten = xlsread('K1.xlsx'); %Eingabe der Messdaten via Excel-Sheet: Spalte 1-3 Eulerwinkel, Spalte 4 x-Pos, Spalte 5 y-Pos
KS_Varianten = xlsread('Varianten_bereinigt.xlsx');

phi1_mess = Messdaten(1:length(Messdaten),1);
theta_mess = Messdaten(1:length(Messdaten),2);
phi2_mess = Messdaten(1:length(Messdaten),3);

phi1_KS = KS_Varianten(1:length(KS_Varianten),1);
theta_KS = KS_Varianten(1:length(KS_Varianten),2);
phi2_KS = KS_Varianten(1:length(KS_Varianten),3);

aa=0;
bb=0;
cc=0;
dd=0;
ee=0;
ff=0;
gg=0;
hh=0;

aufgefuellt = zeros(10,3);
aufgefuellt2 = zeros(10,3);
aufgefuellt3 = zeros(10,3);

aufgefuellt4 = zeros(10,3);
aufgefuellt5 = zeros(10,3);
aufgefuellt6 = zeros(10,3);

Transformation_V(:,:,1)=[0.742	-0.667	-0.075  %alle KS-Varianten nach Kitahara Acta mat 20016
                        0.650	0.742	-0.167
                        0.167	0.075	0.983];
Transformation_V(:,:,2)=[0.075	0.667	-0.742
                        -0.167	0.742	0.650
                        0.983	0.075	0.167];
Transformation_V(:,:,3)=[-0.667	-0.075	0.742
                        0.742	-0.167	0.650
                        0.075	0.983	0.167];
Transformation_V(:,:,4)=[0.667	-0.742	0.075
                        0.742	0.650	-0.167
                        0.075	0.167	0.983];
Transformation_V(:,:,5)=[-0.075	0.742	-0.667
                        -0.167	0.650	0.742
                        0.983	0.167	0.075];
Transformation_V(:,:,6)=[-0.742	0.075	0.667
                        0.650	-0.167	0.742
                        0.167	0.983	0.075];
Transformation_V(:,:,7)=[-0.075	0.667	0.742
                        -0.167	-0.742	0.650
                        0.983	-0.075	0.167];
Transformation_V(:,:,8)=[-0.742	-0.667	0.075
                        0.650	-0.742	-0.167
                        0.167	-0.075	0.983];
Transformation_V(:,:,9)=[0.742	0.075	-0.667
                        0.650	0.167	0.742
                        0.167	-0.983	0.075];
Transformation_V(:,:,10)=[0.075	0.742	0.667
                        -0.167	-0.650	0.742
                        0.983	-0.167	0.075];
Transformation_V(:,:,11)=[-0.667	-0.742	-0.075
                        0.742	-0.650	-0.167
                        0.075	-0.167	0.983];
Transformation_V(:,:,12)=[0.667	-0.075	-0.742
                        0.742	0.167	0.650
                        0.075	-0.983	0.167];
Transformation_V(:,:,13)=[0.667	0.742	-0.075
                        -0.742	0.650	-0.167
                        -0.075	0.167	0.983];
Transformation_V(:,:,14)=[-0.667	0.075	-0.742
                        -0.742	-0.167	0.650
                        -0.075	0.983	0.167];
Transformation_V(:,:,15)=[0.075	-0.667	0.742
                         0.167	0.742	0.650
                        -0.983	0.075	0.167];
Transformation_V(:,:,16)=[0.742	0.667	0.075
                        -0.650	0.742	-0.167
                        -0.167	0.075	0.983];
Transformation_V(:,:,17)=[-0.742	-0.075	-0.667
                        -0.650	-0.167	0.742
                        -0.167	0.983	0.075];
Transformation_V(:,:,18)=[-0.075	-0.742	0.667
                         0.167	0.650	0.742
                        -0.983	0.167	0.075];
Transformation_V(:,:,19)=[0.742	-0.075	0.667
                         0.650	-0.167	-0.742
                         0.167	0.983	-0.075];
Transformation_V(:,:,20)=[0.075	-0.742	-0.667
                        -0.167	0.650	-0.742
                         0.983	0.167	-0.075];
Transformation_V(:,:,21)=[-0.667	0.742	0.075
                         0.742	0.650	0.167
                         0.075	0.167	-0.983];
Transformation_V(:,:,22)=[0.667	0.075	0.742
                         0.742	-0.167	-0.650
                         0.075	0.983	-0.167];
Transformation_V(:,:,23)=[-0.075	-0.667	-0.742
                        -0.167	0.742	-0.650
                         0.983	0.075	-0.167];
Transformation_V(:,:,24)=[-0.742	0.667	-0.075
                         0.650	0.742	0.167
                         0.167	0.075	-0.983];

Transformation_V_inv=[0 0 0];

for i=1:24
        Transformation_V_invH=inv(Transformation_V(:,:,i));
        Transformation_V_inv(1,:,i) = Transformation_V_invH(1,:);
        Transformation_V_inv(2,:,i) = Transformation_V_invH(2,:);
        Transformation_V_inv(3,:,i) = Transformation_V_invH(3,:);
end
                     
Vergleichswert1=[0 0];
Min_gesamt = [0 0];

%%Messdatenauftragung

for j=1:length(Messdaten); 
 
 phi1 = phi1_mess(j);
 theta = theta_mess(j);
 phi2 = phi2_mess(j);
 
 x=phi1;
 y=phi2;
 z=theta;

 Exp=[cos(y)*cos(x)-sin(x)*sin(y)*cos(z) -cos(y)*sin(x)-sin(y)*cos(z)*cos(x) sin(y)*sin(z)
     sin(y)*cos(x)+cos(y)*cos(z)*sin(x) -sin(y)*sin(x)+cos(y)*cos(z)*cos(x) -cos(y)*sin(z)
     sin(z)*sin(x) sin(z)*cos(x) cos(z)];
 
 %Vergleich Reale Messwerte mit idealen KS-Varianten
 for n=1:length(KS_Varianten)
     %Verschiedene Möglichkeiten zum Vergleich
%     winkel1 = acosd(dot(Exp(:,1),Transformation_V(:,1,n)));
%     winkel2 = acosd(dot(Exp(:,2),Transformation_V(:,2,n)));
%     winkel3 = acosd(dot(Exp(:,3),Transformation_V(:,3,n)));
% 
%     Vergleichswert1(j,n) = real(winkel1+winkel2+winkel3);  %Vergleichen der Winkel wenn der minimal wird sollte es passen 
%     %Vergleichswert1(j,n) = real(winkel1)+real(winkel2)+real(winkel3);
    
%     zahl1=dot(Exp(:,1),Transformation_V(:,1,n));
%     zahl2=dot(Exp(:,2),Transformation_V(:,2,n));
%     zahl3=dot(Exp(:,3),Transformation_V(:,3,n));

    Vergleichswert1H=dot(Exp,Transformation_V_inv(:,:,n));
    Vergleichswert1(j,n) = sum(Vergleichswert1H);

%     Vergleichswert1(j,n)=abs(3-(zahl1+zahl2+zahl3)); %Vergleichen der Skalarprodukte (1 ist optimal) wenn der minimal wird sollte es passen
 
 end
 Varianten(j,1) = Messdaten(j,4); %x-Pos für Bild
 Varianten(j,2) = Messdaten(j,5); %y-Pos für Bild
 [Min,Pos] = min(Vergleichswert1(j,:));
 Min_gesamt(j) = Min;
 
 if Min<1.9    
     Varianten(j,3) = Pos;
 else
     Varianten(j,3) = 0;
 end
 
%Stereographische Projektion 
for i = 1:size(v,1);
     if v(i,3) < 0;
         v = v(i,:) * (-1);     
     end
     
     rot_basis = [v(i,:); 0 0 0; 0 0 0]';
     rot_basis = eulerAngleRotation( phi1, theta, phi2, rot_basis );
     v_rot(i,:) = rot_basis(:,1);
     v_rot(i,:) = v_rot(i,:)/ norm( v_rot(i,:) );

     if v_rot(i,3) < 0;
         v_rot(i,:) = v_rot(i,:) * (-1);
     end
     
     t = 1 / (v_rot(i,3) +1);
     v_proj_temp(i,:) = [ v_rot(i,1) * t, v_rot(i,2)*t ];

end

%Zeilenweise zusammenstellung des v_proj - Vektors aus allen Messdaten

v_proj(1,:,j) = v_proj_temp(1,:);
v_proj(2,:,j) = v_proj_temp(2,:);
v_proj(3,:,j) = v_proj_temp(3,:);

end

%%KS-Variantenauftragung

for j=1:length(KS_Varianten);
 
 phi1 = phi1_KS(j);
 theta = theta_KS(j);
 phi2 = phi2_KS(j);
    
for i = 1:size(v,1);
     if v(i,3) < 0;
         v = v(i,:) * (-1); 
     end
     
     rot_basis = [v(i,:); 0 0 0; 0 0 0]';
     rot_basis = eulerAngleRotation( phi1, theta, phi2, rot_basis );
     v_rot(i,:) = rot_basis(:,1);
     v_rot(i,:) = v_rot(i,:)/ norm( v_rot(i,:) );
     
     if v_rot(i,3) < 0;
         v_rot(i,:) = v_rot(i,:) * (-1);
     end
     t = 1 / (v_rot(i,3) +1); 
     v_proj_temp(i,:) = [ v_rot(i,1) * t, v_rot(i,2)*t ];

end

v_proj_KS(1,:,j) = v_proj_temp(1,:);
v_proj_KS(2,:,j) = v_proj_temp(2,:);
v_proj_KS(3,:,j) = v_proj_temp(3,:);

end

%Herstellung der richtigen Größe des Arrays für den nächsten Schritt
for o=1:length(Messdaten);
    for p=1:length(v);
        v_proj(p,3,o)=0;
    end    
end

%Darstelung der Messdaten in der Polfigur
myPolar([0 0],[-1,1],1, ['k', '.']);
hold all
myPolar([pi/2 pi/2],[-1,1],1, ['k', '.']);
hold all
title('Vergleich der Messdaten mit idealen KS-Varianten')

for k=1:length(Messdaten);
plot( v_proj(:,1,k), v_proj(:,2,k), '.', 'MarkerSize', 1, 'MarkerFaceColor','b'); %v_rot(:,1), v_rot(:,2), '*' ); %v(:,1), v(:,2), '*'); %
axis([-1 1 -1 1]);
view(270,90);
hold on;
end

%Überlagerung mit dem idealen KS-Pattern
for l=1:length(KS_Varianten);
plot( v_proj_KS(:,1,l), v_proj_KS(:,2,l), 'o', 'MarkerSize', 5, 'MarkerFaceColor','none','MarkerEdgeColor','k'); %v_rot(:,1), v_rot(:,2), '*' ); %v(:,1), v(:,2), '*'); %
axis([-1 1 -1 1]);
view(270,90);
hold on;
labels = cellstr( num2str([l]') );  
text(v_proj_KS(:,1,l), v_proj_KS(:,2,l),labels);
end        

%Umwandeln der x-,y-Koordinaten der Varianten in darstellbares Bildformat
y=1;
x=1;
for s=1:length(Varianten)-1;   
    if Varianten(s,2) == Varianten(s+1,2);
        Bild(y,x) = Varianten(s,3);
        x=x+1;    
    else
        Bild(y,x) = Varianten(s,3);
        y=y+1;
        x=1;
    end
end
Bild_original = Bild+1; %"Variante 0" (ungültige Messpunkte) erhalten Zahlenwert 1 für die Ausgabe als Bild
Groesse_Bild = size(Bild);


%eigene Colorbar für Varianten auswählen
 
  E_Color = [1.00 1.00 1.00; %D: Varianten des ersten Paketes eingefärbt
     1.00 0.00 0.00; %V1
     0.00 1.00 0.00; 
     0.00 0.00 1.00;
     0.70 0.00 0.00; 
     0.00 0.70 0.00; 
     0.00 0.00 0.70; 
     
     0.55 0.55 0.55; %V7 
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77; 
     
     0.55 0.55 0.55; %V13
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77; 
     
     0.55 0.55 0.55; %V19
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77]; 
 
F_Color = [1.00 1.00 1.00; %D: Varianten des ersten Paketes eingefärbt
    
     0.55 0.55 0.55; %V1
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     
     1.00 0.00 0.00; %V7
     0.00 1.00 0.00; 
     0.00 0.00 1.00;
     0.70 0.00 0.00; 
     0.00 0.70 0.00; 
     0.00 0.00 0.70;
     
     0.55 0.55 0.55; %V13     
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77; 
     
     0.55 0.55 0.55; %V19
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77];  
 
 G_Color = [1.00 1.00 1.00; %D: Varianten des ersten Paketes eingefärbt
    
     0.55 0.55 0.55; %V1
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     
     0.55 0.55 0.55; %V7
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     
     1.00 0.00 0.00; %V13
     0.00 1.00 0.00; 
     0.00 0.00 1.00;
     0.70 0.00 0.00; 
     0.00 0.70 0.00; 
     0.00 0.00 0.70; 
     
     0.55 0.55 0.55; %V19
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77];  
 
 H_Color = [1.00 1.00 1.00; %D: Varianten des ersten Paketes eingefärbt
    
     0.55 0.55 0.55; %V1
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     
     0.55 0.55 0.55; %V7
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     
     0.55 0.55 0.55; %V13
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     0.55 0.55 0.55;
     0.66 0.66 0.66; 
     0.77 0.77 0.77;
     
     1.00 0.00 0.00; %V19
     0.00 1.00 0.00; 
     0.00 0.00 1.00;
     0.70 0.00 0.00; 
     0.00 0.70 0.00; 
     0.00 0.00 0.70];  

%Darstellung des Bildes 
imtool(Bild+1,E_Color)
title('Varianten innerhalb eines Paketes 1')
imtool(Bild+1, F_Color)
title('Varianten innerhalb eines Paketes 2')
imtool(Bild+1, G_Color)
title('Varianten innerhalb eines Paketes 3')
imtool(Bild+1, H_Color)
title('Varianten innerhalb eines Paketes 4')
