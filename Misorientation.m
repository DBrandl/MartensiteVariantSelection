V = [ 0.742	-0.667	-0.075 %verschiebungsvektor der Variante
0.65	0.742	-0.167
0.167	0.075	0.983 ];

a=1; %verkippungswinkel in grad
n1=-1/sqrt(2); %verkippungsachse als einheitsvektor!!!!
n2=0/sqrt(2);
n3=1/sqrt(2);

g11=n1^2*(1-cosd(a))+cosd(a);
g12=n1*n2*(1-cosd(a))-n3*sind(a);
g13=n1*n3*(1-cosd(a))+n2*sind(a);
g21=n2*n1*(1-cosd(a))+n3*sind(a);
g22=n2^2*(1-cosd(a))+cosd(a);
g23=n2*n3*(1-cosd(a))-n1*sind(a);
g31=n3*n1*(1-cosd(a))-n2*sind(a);
g32=n3*n2*(1-cosd(a))+n1*sind(a);
g33=n3^2*(1-cosd(a))+cosd(a);

G=[g11 g12 g23
    g21 g22 g23
    g31 g32 g33];

T=transpose(G);

M=V*G;

PHI=acos(M(3,3)); %in Radianten

if M(3,2) < 0
    phi1=pi-asin(M(3,1)/sin(PHI));
else 
    phi1=-(sign(M(3,1))-sign(M(3,2)))*pi+asin(M(3,1)/sin(PHI));
end

if M(2,3) < 0
    phi2=-(sign(M(1,3))+sign(M(2,3)))*pi+asin(M(1,3)/sin(PHI));
else 
    phi2=pi-asin(M(1,3)/sin(PHI));
end


