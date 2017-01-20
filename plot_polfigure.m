v = [1 1 0; % var 4
     1 0 1; % var 1
     0 1 1; % var 5
     -1 1 0; % var 2 
     -1 0 1; % var 3
     0 -1 1]; % var 6
 
 v_rot = [0 0 0];
 v_proj = [0 0];
 
 phi1 = (pi/180) * 182;
 theta = (pi/180) * 74;
 phi2 = (pi/180) * 306;
    
for i = 1:size(v,1)
     if v(i,3) < 0
         v = v(i,:) * (-1);
     end
     
     rot_basis = [v(i,:); 0 0 0; 0 0 0]';
     rot_basis = eulerAngleRotation( phi1, theta, phi2, rot_basis );
     v_rot(i,:) = rot_basis(:,1);
     v_rot(i,:) = v_rot(i,:)/ norm( v_rot(i,:) );
     
     % For the pole figure only one half of the sphere is important,
     % therefore if the vector cuts the unit sphere in the lower half with
     % z<0 then the direction of the vector is inversed.
     % (coordinates are adopted according to the point symmetry around the point [0 0 0].
     if v_rot(i,3) < 0
         v_rot(i,:) = v_rot(i,:) * (-1);
     end
     
     % Jetzt wird die Stereographische Projection ermittelt. Dazu wird zuerst
     % der Schnitt mit der Einheitskugel berechnet. Die Stereographische
     % Projection ist nun der Schnittpunkt der Gerade definiert durch den Punkt auf der Einheitskugel
     % und dem Punkt O [0 0 -1] mit der Ebene z=0.
     t = 1 / (v_rot(i,3) +1); % diese Gleichung kommt aus der Parameterform für die Z-Komponente
     v_proj(i,:) = [ v_rot(i,1) * t, v_rot(i,2)*t ];

end


myPolar([0 0],[-1,1],1, ['k', '.']);
hold all
myPolar([pi/2 pi/2],[-1,1],1, ['k', '.']);
hold all
plot( v_proj(:,1), v_proj(:,2), 'o', 'MarkerSize', 10, 'MarkerFaceColor','b'); %v_rot(:,1), v_rot(:,2), '*' ); %v(:,1), v(:,2), '*'); %
axis([-1 1 -1 1]);
% angles = linspace(0,2*pi);
% plot( cos(angles), sin(angles),'b');
view(90,90);


