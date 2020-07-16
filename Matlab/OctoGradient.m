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
A = (1 + gradAG);
G = (1 - gradAG);
D = (1 + gradDF);
F = (1 - gradDF);
C = (1 + gradCE);
E = (1 - gradCE);
B = (1 + gradBH);
H = (1 - gradBH);

% Select randon number between 0 and 8 (8 being the total chance given the
% fact that each face had chance of 1 to start. This was mostly chosen so 
% the math requires less operations (no need to multiply the scaled 
% gradients with 0.8))
randomNr = rand * 8;
face_probability = -1;

% Prepare two random values to chose random point on triangle(=face)
S = rand;
T = rand;
if S + T > 1
    S = 1 - S;
    T = 1 - T;
end

% Determine face and create direction vector
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
        
% calculate probability
% Pdf = face_probability * 1/8 (=uniform chance of one face) / area of face (=4*PI/8) -> reduces to face_probability * INV_$_PI_F
Pdf = face_probability * (1 / (4*pi));
F = Kd * (1 / (4*pi));
end

