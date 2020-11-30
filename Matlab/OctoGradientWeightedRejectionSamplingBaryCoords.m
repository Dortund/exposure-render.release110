function [p_F, Wi, Pdf, Face, Tries] = OctoGradientWeightedRejectionSamplingBaryCoords(Kd, datapoints)
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
%p = 0.00387597;
f_AEHD = f_AEHD^p;
f_GCBF = f_GCBF^p;
f_ABFE = f_ABFE^p;
f_GHDC = f_GHDC^p;
f_ABCD = f_ABCD^p;
f_GHEF = f_GHEF^p;

%sum = AEHD + GCBF + ABFE +GHDC + ABCD+GHEF
%{
faces = [
    [f_GHEF, f_GCBF, f_GHDC]; ... %G
    [f_GHEF, f_ABFE, f_GCBF]; ... %F
    [f_GHEF, f_AEHD, f_ABFE]; ... %E
    [f_GHEF, f_GHDC, f_AEHD]; ... %H
    [f_ABCD, f_GCBF, f_GHDC]; ... %C
    [f_ABCD, f_ABFE, f_GCBF]; ... %B
    [f_ABCD, f_AEHD, f_ABFE]; ... %A
    [f_ABCD, f_GHDC, f_AEHD]      %D
    ];
%}
faces = [
    [f_GCBF, f_GHDC, f_GHEF]; ... %G
    [f_GCBF, f_GHEF, f_ABFE]; ... %F
    [f_GCBF, f_ABFE, f_ABCD]; ... %E
    [f_GCBF, f_ABCD, f_GHDC]; ... %H
    [f_ABCD, f_GHDC, f_GHEF]; ... %C
    [f_ABCD, f_GHEF, f_ABFE]; ... %B
    [f_ABCD, f_ABFE, f_ABCD]; ... %A
    [f_ABCD, f_ABCD, f_GHDC]      %D
    ];
sum = 0;
for i = 1:8
    sumFace = (1/2)*((pi-2)*faces(i,1)+faces(i,2)+faces(i,3));
    sum = sum + sumFace;
end


accepted = false;
Tries = 0;

while (not(accepted))
    Tries = Tries + 1;
%     temp = rand();
%     if (temp <= 1/3)
%         rx = 0;
%     elseif temp <= 2/3
%         rx = 0.5;
%     else
%         rx = 1;
%     end
%     
%     temp = rand();
%     if temp <= 0.25
%         ry = 0;%1/8;
%     elseif temp <= 0.5
%         ry = 1/4;%1/8+1/4;
%     elseif temp <= 0.75
%         ry = 1/2;%1/8+2/4;
%     else
%         ry = 3/4;%1/8+3/4;
%     end
    
    rx = rand();
    ry = rand();

    z = 1 - 2 * rx;
    r = sqrt(max(0, 1 - z*z));
    phi = 2 * pi * ry;
    x = r * cos(phi);
    y = r * sin(phi);
    
    if (z > 0)
        if (phi <= pi/2)    
%             A = GHEF;
%             C = GCBF;
%             B = GHDC;
                p_A = f_GCBF;
					p_B = f_GHDC;
					p_C = f_GHEF;
            Face = [0,0,0,0,0,0,1,0];
        elseif (phi <= pi)
%             A = GHEF;
%             C = ABFE;
%             B = GCBF;
p_A = f_GCBF;
					p_B = f_GHEF;
					p_C = f_ABFE;
            Face = [0,0,0,0,0,1,0,0];
        elseif (phi <= 3*pi/2)
%             A = GHEF;
%             C = AEHD;
%             B = ABFE;
p_A = f_GCBF;
					p_B = f_ABFE;
					p_C = f_ABCD;
            Face = [0,0,0,0,1,0,0,0];
        else
%             A = GHEF;
%             C = GHDC;
%             B = AEHD;
p_A = f_GCBF;
					p_B = f_ABCD;
					p_C = f_GHDC;
            Face = [0,0,0,0,0,0,0,1];
        end
    else
        if (phi <= pi/2)    
%             A = ABCD;
%             C = GCBF;
%             B = GHDC;
p_A = f_AEHD;
					p_B = f_GHDC;
					p_C = f_GHEF;
            Face = [0,0,1,0,0,0,0,0];
        elseif (phi <= pi)
%             A = ABCD;
%             C = ABFE;
%             B = GCBF;
p_A = f_AEHD;
					p_B = f_GHEF;
					p_C = f_ABFE;
            Face = [0,1,0,0,0,0,0,0];
        elseif (phi <= 3*pi/2)
%             A = ABCD;
%             C = AEHD;
%             B = ABFE;
p_A = f_AEHD;
					p_B = f_ABFE;
					p_C = f_ABCD;
            Face = [1,0,0,0,0,0,0,0];
        else
%             A = ABCD;
%             C = GHDC;
%             B = AEHD;
p_A = f_AEHD;
					p_B = f_ABCD;
					p_C = f_GHDC;
            Face = [0,0,0,1,0,0,0,0];
        end
    end
    
    % random uniform direction: x
    % take dot product for x and every corner vector
    % clamp those in range [0,1]
    % multiply by corner weights, sum these values for PDF chance
      
    %Pdf = (((C-B)*2*rem(phi, pi/2)/pi) -A)*2*(rx*pi/2)/pi+A;
    %Pdf = (((C-B)*rem(ry, 1/4)*4 + B) - A)*rem(rx, 1/2)*2 + A;
    party = rem(ry, 1/4)*4;
    partx = 1-abs(rx-(1/2))*2;%rem(rx, 1/2)*2;
    hor = (p_C-p_B)*party+p_B;
    ver = (hor-p_A)*partx+p_A;
    Pdf = ver/sum;
    %{
    vec = [x, y ,z];
    
    aehd = min(max(dot([0, 0, -1], vec), 0), 1) * f_AEHD;
	gcbf = min(max(dot([0, 0,  1], vec), 0), 1) * f_GCBF;
	ghef = min(max(dot([0,  1, 0], vec), 0), 1) * f_GHEF;
	abcd = min(max(dot([0, -1, 0], vec), 0), 1) * f_ABCD;
	abfe = min(max(dot([-1, 0, 0], vec), 0), 1) * f_ABFE;
	ghdc = min(max(dot([ 1, 0, 0], vec), 0), 1) * f_GHDC;
    
    p = aehd + gcbf + ghef + abcd + abfe + ghdc;
	Pdf2 = p / (abs(x)+abs(y)+abs(z)) / sum;
    %Pdf = Pdf2;
    %}
    %accepted = rand() <= Pdf;
    %accepted = rand() <= Pdf*Pdf;
    
    maxP = max(f_AEHD, max(f_GCBF, max(f_ABFE, max(f_GHDC, max(f_ABCD, f_GHEF)))));
    accepted = rand() <= ver / maxP;

end

%Wi = [x, y, z] * Pdf2;
Wi = [x, y, z] * Pdf;
    

end


function [Wi, Pdf] = WeightedTriangle(a, b, c)
% VertexValues = [A, B, C]
%
%       C
%      /\_ _ D
%     /  \
%   A/____\B

%a = VertexValues(1);
%b = VertexValues(2);
%c = VertexValues(3);

cMinB = c-b;
A = 0.5*(b+c);

% y = cofA * x^2 + cofB * x + cofC
cofA = cMinB / 2;
cofB = b;
cofC = 0;

bcVector = [-0.5, sqrt(3)/2];

bc = -1;
if b == c
    bc = rand;
else
    bc = (-cofB + sqrt(cofB.^2 - 4*cofA*(cofC - (rand*A)))) / (2*cofA);
end

d = (cMinB*bc + b);

adVector = (bcVector * bc) + [1, 0];


dMinA = d-a;
A = 0.5*(a+d);

cofA2 = dMinA / 2;
cofB2 = a;
cofC2 = 0;

ad = -1;
pU = sqrt(rand);

if a == d
    ad = pU;
else
    ad = (-cofB2 + sqrt(cofB2^2 - 4*cofA2*(cofC2 - (pU*A)))) / (2*cofA2);
end


Wi = ad * adVector;
Pdf = (dMinA*ad + a);

%th = -1/12*2*pi; %30 degrees clockwise
th = 2/12*2*pi; %60 degrees anti-clockwise
rm = [[cos(th) -sin(th)];[sin(th) cos(th)]];

Wi = transpose(rm*transpose(Wi));
end