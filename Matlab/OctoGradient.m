function [F, Wi, Pdf, Face] = OctoGradient(Kd, datapoints)
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

gradAG = datapoints(7) - datapoints(1);
gradDF = datapoints(6) - datapoints(4);
gradCE = datapoints(5) - datapoints(3);
gradBH = datapoints(8) - datapoints(2);

gradAG = gradAG / 101 * .99;
gradDF = gradDF / 101 * .99;
gradCE = gradCE / 101 * .99;
gradBH = gradBH / 101 * .99;

A = (1 + gradAG);
G = (1 - gradAG);
D = (1 + gradDF);
F = (1 - gradDF);
C = (1 + gradCE);
E = (1 - gradCE);
B = (1 + gradBH);
H = (1 - gradBH);

randomNr = rand * 8;
face_probability = -1;

S = rand;
T = rand;
if S + T > 1
    S = 1 - S;
    T = 1 - T;
end

if randomNr < E
    v0 = [0, 0, -1];
    v1 = [0, 1, 0] - v0;
    v2 = [-1, 0, 0] - v0;
    Wi = v0 + S * v1 + T * v2;
    Wi = Wi / norm(Wi);
    face_probability = E;
    Face = [0,0,0,0,1,0,0,0];
elseif randomNr < E + H
    v0 = [1, 0, 0];
    v1 = [0, 1, 0] - v0;
    v2 = [0, 0, -1] - v0;
    Wi = v0 + S * v1 + T * v2;
    Wi = Wi / norm(Wi);
    face_probability = H;
    Face = [0,0,0,0,0,0,0,1];
elseif randomNr < E + H + G
    v0 = [0, 0, 1];
    v1 = [0, 1, 0] - v0;
    v2 = [1, 0, 0] - v0;
    Wi = v0 + S * v1 + T * v2;
    Wi = Wi / norm(Wi);
    face_probability = G;
    Face = [0,0,0,0,0,0,1,0];
elseif randomNr < E + H + G + F
    v0 = [-1, 0, 0];
    v1 = [0, 1, 0] - v0;
    v2 = [0, 0, 1] - v0;
    Wi = v0 + S * v1 + T * v2;
    Wi = Wi / norm(Wi);
    face_probability = F;
    Face = [0,0,0,0,0,1,0,0];
elseif randomNr < E + H + G + F + A
    v0 = [0, 0, -1];
    v1 = [0, -1, 0] - v0;
    v2 = [-1, 0, 0] - v0;
    Wi = v0 + S * v1 + T * v2;
    Wi = Wi / norm(Wi);
    face_probability = A;
    Face = [1,0,0,0,0,0,0,0];
elseif randomNr < E + H + G + F + A + D
    v0 = [1, 0, 0];
    v1 = [0, -1, 0] - v0;
    v2 = [0, 0, -1] - v0;
    Wi = v0 + S * v1 + T * v2;
    Wi = Wi / norm(Wi);
    face_probability = D;
    Face = [0,0,0,1,0,0,0,0];
elseif randomNr < E + H + G + F + A + D + C
    v0 = [0, 0, 1];
    v1 = [0, -1, 0] - v0;
    v2 = [1, 0, 0] - v0;
    Wi = v0 + S * v1 + T * v2;
    Wi = Wi / norm(Wi);
    face_probability = C;
    Face = [0,0,1,0,0,0,0,0];
else
    v0 = [-1, 0, 0];
    v1 = [0, -1, 0] - v0;
    v2 = [0, 0, 1] - v0;
    Wi = v0 + S * v1 + T * v2;
    Wi = Wi / norm(Wi);
    face_probability = B;
    Face = [0,1,0,0,0,0,0,0];
end
        
Pdf = face_probability * (1 / (4*pi));
F = Kd * (1 / (4*pi));
end

