function [p_F, Wi, Pdf, Face, Phi, Theta] = OctoGradientWeightedInverseSampling(Kd, datapoints)
%OCTOGRADIENT Summary of this function goes here
%   Detailed explanation goes here
%{
    [I]=G________F  +y
	    /      /|    ^
	   /      / |    |
	 H/_____E/  |
	  |  C  |  /B       +z
	  |     | /        /
	  |_____|/       -z
	 D     A=(0,0,0)
	    +x<->-x
%}

% Calculate gradients between opposing faces
gradAG = datapoints(7) - datapoints(1);
gradDF = datapoints(6) - datapoints(4);
gradCE = datapoints(5) - datapoints(3);
gradBH = datapoints(8) - datapoints(2);

% Scale gradients between -1 and and 1 using know max gradient (101)
gradAG = gradAG / 101;
gradDF = gradDF / 101;
gradCE = gradCE / 101;
gradBH = gradBH / 101;

% Offset uniform chance of choosing a face
p_A = (1 + gradAG);
p_G = (1 - gradAG);
p_D = (1 + gradDF);
p_F = (1 - gradDF);
p_C = (1 + gradCE);
p_E = (1 - gradCE);
p_B = (1 + gradBH);
p_H = (1 - gradBH);

f_AEHD = (p_A+p_E+p_H+p_D) / 4;
f_GCBF = (p_G+p_C+p_B+p_F) / 4;
f_ABFE = (p_A+p_B+p_F+p_E) / 4;
f_GHDC = (p_G+p_H+p_D+p_C) / 4;
f_ABCD = (p_A+p_B+p_C+p_D) / 4;
f_GHEF = (p_G+p_H+p_E+p_F) / 4;

p = 1;
f_AEHD = f_AEHD^p;
f_GCBF = f_GCBF^p;
f_ABFE = f_ABFE^p;
f_GHDC = f_GHDC^p;
f_ABCD = f_ABCD^p;
f_GHEF = f_GHEF^p;

faces = [
    [f_GCBF, f_GHDC, f_GHEF]; ... %G
    [f_GCBF, f_GHEF, f_ABFE]; ... %F
    [f_GCBF, f_ABFE, f_ABCD]; ... %E
    [f_GCBF, f_ABCD, f_GHDC]; ... %H
    [f_AEHD, f_GHDC, f_GHEF]; ... %C
    [f_AEHD, f_GHEF, f_ABFE]; ... %B
    [f_AEHD, f_ABFE, f_ABCD]; ... %A
    [f_AEHD, f_ABCD, f_GHDC]      %D
    ];
sum = zeros(8,1);
for i = 1:8
    sumFace = (1/2)*((pi-2)*faces(i,1)+faces(i,2)+faces(i,3));
    if i == 1
        sum(i) = sumFace;
    else
        sum(i) = sum(i-1) + sumFace;
    end
end

faceRnd = rand() * sum(8);
%faceRnd = sum(7);

if faceRnd <= sum(1) %G
   A = faces(1,3);
   B = faces(1,2);
   C = faces(1,1);
   offsetPhi = 1;
   offsetTheta = 0;
   Face = [0,0,0,0,0,0,1,0];
elseif faceRnd <= sum(2) %F
   A = faces(2,3);
   B = faces(2,2);
   C = faces(2,1);
   offsetPhi = 1;
   offsetTheta = pi/2;
   Face = [0,0,0,0,0,1,0,0];
elseif faceRnd <= sum(3) %E
   A = faces(3,3);
   B = faces(3,2);
   C = faces(3,1);
   offsetPhi = 1;
   offsetTheta = pi;
   Face = [0,0,0,0,1,0,0,0];
elseif faceRnd <= sum(4) %H
   A = faces(4,3);
   B = faces(4,2);
   C = faces(4,1);
   offsetPhi = 1;
   offsetTheta = 3*pi/2;
   Face = [0,0,0,0,0,0,0,1];
elseif faceRnd <= sum(5) %C
   A = faces(5,3);
   B = faces(5,2);
   C = faces(5,1);
   offsetPhi = -1;
   offsetTheta = 0;
   Face = [0,0,1,0,0,0,0,0];
elseif faceRnd <= sum(6) %B
   A = faces(6,3);
   B = faces(6,2);
   C = faces(6,1);
   offsetPhi = -1;
   offsetTheta = pi/2;
   Face = [0,1,0,0,0,0,0,0];
elseif faceRnd <= sum(7) %A
   A = faces(7,3);
   B = faces(7,2);
   C = faces(7,1);
   offsetPhi = -1;
   offsetTheta = pi;
   Face = [1,0,0,0,0,0,0,0];
elseif faceRnd <= sum(8) %D
   A = faces(8,3);
   B = faces(8,2);
   C = faces(8,1);
   offsetPhi = -1;
   offsetTheta = 3*pi/2;
   Face = [0,0,0,1,0,0,0,0];
end

%A = 5;
%B = 10;
%C = 2;

xrnd = rand();
yrnd = rand();

%xrnd = 1;
%yrnd = 0;

sumFace = (1/2)*((pi-2)*C+A+B);

% inverse integrate (B-A)x+A dx, x=0..x
if (A ~= B)
    xScaled = xrnd * sumFace;
    %xInv = (-2*pi*B - pi^2*C + 2*pi*C)/(4*(A - B)) + (pi*sqrt(8*A*xScaled + 4*B^2 - 8*B*C + 4*pi*B*C - 8*B*xScaled + 4*C^2 + pi^2*C^2 - 4*pi*C^2))/(4*(A - B));
    xInv = (pi * (sqrt(8 * A*xScaled + (2*B + (pi - 2)*C)^2 - 8*B*xScaled) - 2*B - (pi - 2)*C)) / (4*(A - B));
else
    xInv = xrnd * (pi/2);
end

D = (A-B)*(xInv/(pi/2))+B;

if (C ~= D)
    yScaled = yrnd * (C+D)/2;
    yInv = (C - sqrt(C^2 - 2*C*yScaled + 2*D*yScaled))/(C - D);
else
    yInv = yrnd;
end

theta = xInv + offsetTheta;
phi = acos(1 - yInv);
x = sin(phi) * cos(theta);
y = sin(phi) * sin(theta);
z = cos(phi) * offsetPhi;
Wi = [x, y, z];

yPart = phi / (0.5*pi);
P = (D-C)*yPart+C;
Pdf = P/sum(8);

check = CheckPdf(Wi, datapoints);
%if Pdf ~= check
if abs(Pdf-check) > 0.00001
    asdf = 9
end

if offsetPhi == 1
    Phi = phi;
else
    Phi = pi - phi;
end
Theta = theta;
    
%if x < 0.2 && y < 0.2 && z < 0.2
    fin = -9;
%end

end

function Pdf = CheckPdf(Wi, datapoints)
   % Calculate gradients between opposing faces
    gradAG = datapoints(7) - datapoints(1);
    gradDF = datapoints(6) - datapoints(4);
    gradCE = datapoints(5) - datapoints(3);
    gradBH = datapoints(8) - datapoints(2);

    % Scale gradients between -1 and and 1 using know max gradient (101)
    gradAG = gradAG / 101;
    gradDF = gradDF / 101;
    gradCE = gradCE / 101;
    gradBH = gradBH / 101;

    % Offset uniform chance of choosing a face
    p_A = (1 + gradAG);
    p_G = (1 - gradAG);
    p_D = (1 + gradDF);
    p_F = (1 - gradDF);
    p_C = (1 + gradCE);
    p_E = (1 - gradCE);
    p_B = (1 + gradBH);
    p_H = (1 - gradBH);

    f_AEHD = (p_A+p_E+p_H+p_D) / 4;
    f_GCBF = (p_G+p_C+p_B+p_F) / 4;
    f_ABFE = (p_A+p_B+p_F+p_E) / 4;
    f_GHDC = (p_G+p_H+p_D+p_C) / 4;
    f_ABCD = (p_A+p_B+p_C+p_D) / 4;
    f_GHEF = (p_G+p_H+p_E+p_F) / 4;

    p = 1;
    f_AEHD = f_AEHD^p;
    f_GCBF = f_GCBF^p;
    f_ABFE = f_ABFE^p;
    f_GHDC = f_GHDC^p;
    f_ABCD = f_ABCD^p;
    f_GHEF = f_GHEF^p;

    faces = [
        [f_GCBF, f_GHDC, f_GHEF]; ... %G
        [f_GCBF, f_GHEF, f_ABFE]; ... %F
        [f_GCBF, f_ABFE, f_ABCD]; ... %B
        [f_GCBF, f_ABCD, f_GHDC]; ... %C
        [f_AEHD, f_GHDC, f_GHEF]; ... %H
        [f_AEHD, f_GHEF, f_ABFE]; ... %E
        [f_AEHD, f_ABFE, f_ABCD]; ... %A
        [f_AEHD, f_ABCD, f_GHDC]      %D
        ];
    
    sum = 0;
    for i = 1:8
        sum = sum + (1/2)*((pi-2)*faces(i,1)+faces(i,2)+faces(i,3));
    end
    
    theta = atan(Wi(2) / Wi(1));
    phi = acos(Wi(3));

    if Wi(3) > 0 
        if (theta <= pi / 2)
            face = 0;
        elseif (theta <= pi)
            face = 1;
        elseif (theta <= 3 * pi / 2)
            face = 2;
        else
            face = 3;
        end
    else
        if (theta <= pi / 2)
            face = 4;
        elseif (theta <= pi)
            face = 5;
            elseif (theta <= 3 * pi / 2)
            face = 6;
        else
            face = 7;
        end
    end

    t_A = faces(face+1,3);
    t_B = faces(face+1,2);
    t_C = faces(face+1,1);

    rx = theta / (2*pi);
    ry = phi / pi;

    partx = (rx - mod(face, 4) * 0.25) * 4;
    party = 1 - abs(ry - 0.5) * 2;
    hor = (t_A - t_B)*partx + t_B;
    ver = (hor - t_C)*party + t_C;
    Pdf = ver / sum;
end