Vbase = [1/sqrt(3),1/sqrt(3),1/sqrt(3)];
%Vtemp = [-1/sqrt(3),1/sqrt(3),1/sqrt(3)];

Vx = [2/sqrt(3),0,0];
Vy = [0,2/sqrt(3),0];
Vz = [0,0,2/sqrt(3)];

n = 1000;
Vrand = zeros(n,3);
Vnorm = zeros(n,3);
Vrand2 = zeros(n,3);
Vnorm2 = zeros(n,3);
Vrand3 = zeros(n,3);
Vnorm3 = zeros(n,3);

Vreal = zeros(n,3);

for i = 1:n
    Vnew = Vbase - Vx*rand - Vy*rand;
    Vrand(i,:) = Vnew;
    Vnorm(i,:) = Vnew / norm(Vnew);
    
    Vnew = Vbase - Vy*rand - Vz*rand;
    Vrand2(i,:) = Vnew;
    Vnorm2(i,:) = Vnew / norm(Vnew);
    
    Vnew = Vbase - Vx*rand - Vz*rand;
    Vrand3(i,:) = Vnew;
    Vnorm3(i,:) = Vnew / norm(Vnew);
    
    %z = (1 - 2 * rand);
    %z = ((1/sqrt(3)) - (2/sqrt(3)) * rand);
	%r = sqrt(max(0, 1 - z*z));
	%phi = 2 * pi * rand / 4 + 0.25*pi;
	%x = r * cos(phi);
	%y = r * sin(phi);
	
    theta = 2 * pi * rand / 4 + 0.25*pi;
    phi = acos((1/sqrt(3)) - (2/sqrt(3)) * rand);
    x = sin(phi) * cos(theta);
    y = sin(phi) * sin(theta);
    z = cos(phi);
    
    Vnew = [x, y, z];
    Vreal(i,:) = Vnew;
end

figure;
scatter3(Vrand(:,1),Vrand(:,2),Vrand(:,3),'b');
hold on;
scatter3(Vnorm(:,1),Vnorm(:,2),Vnorm(:,3),'r');

scatter3(Vrand2(:,1),Vrand2(:,2),Vrand2(:,3),'g');
scatter3(Vnorm2(:,1),Vnorm2(:,2),Vnorm2(:,3),'c');

scatter3(Vrand3(:,1),Vrand3(:,2),Vrand3(:,3),'y');
scatter3(Vnorm3(:,1),Vnorm3(:,2),Vnorm3(:,3),'m');

scatter3(Vreal(:,1),Vreal(:,2),Vreal(:,3),'k');
%hold off;