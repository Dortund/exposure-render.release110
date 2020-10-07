n = 1000000;
datapoints = [1,2,2];
dirs = zeros(n,2);

for i = 1:n
    [wi, bc, pdf] = WeightedTriangle(datapoints);
    dirs(i,:) = wi;
end
%{
figure;
scatter(dirs(:,1), dirs(:,2));
xlim([0 1])
ylim([-1 1])
%}

figure;
histogram(dirs(:,1),150,'Normalization','probability');
%histogram(dirs(:,1),'Normalization','probability');

%{
[X, Y] = meshgrid(0:0.1:1,1:1:1000);
Z = (-1 + sqrt(1^2 - 4.*((Y-1) ./ 2) .* (0 - (X .* 0.5.*(1+Y))))) ./ (2 .* ((Y-1) ./ 2));
figure;
surf(X,Y,Z);
%}
%{
[X, Y] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
Z = (-1 + sqrt(1^2 - 4.*((Y-1) ./ 2) .* (0 - (X .* 0.5.*(1+Y))))) ./ (2 .* ((Y-1) ./ 2));
Q = ((((Y-1) ./ 2).*(X.^2) + 1.*X + 0)) ./ (0.5.*(1+Y));
figure;
hold on;
for i = 1:size(Y,1)
    plot(X(i,:),Z(i,:));
    plot(X(i,:),Q(i,:));
end
T = X.^2;
R = sqrt(X);
plot(X(1,:),T,'g--o');
plot(X(1,:),R,'b--o');
hold off;
xlim([0 1]);
ylim([0 1]);
%}
%{
[X3, Y3] = meshgrid(0:0.005:1,1:0.05:2);
N = X3.^(1./Y3);
figure;
hold on;
for i = 1:size(Y3,1)
    plot(X3(i,:),N(i,:));
end
hold off;
xlim([0 1]);
ylim([0 1]);
%}
%{
[X4, Y4] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
%a = 2 ./ (1 + Y4)
%b = Y4 .* a
%((b - a) / 2) * X.^2 + a .* X;
M = (((Y4 .* (2 ./ (1 + Y4))) - (2 ./ (1 + Y4))) ./ 2) .* X4.^2 + (2 ./ (1 + Y4)) .* X4;
figure;
hold on;
for i = 1:size(Y4,1)
    plot(X4(i,:),M(i,:));
end
T = X4.^2;
plot(X4(1,:),T,'g--o');
hold off;
xlim([0 1]);
ylim([0 1]);
%}
%{
[X5, Y5] = meshgrid(0:0.005:1,1:1:100);
B = X5.^Y5;
V = X5.^(1./Y5);
figure;
hold on;
for i = 1:size(Y5,1)
    plot(X5(i,:),B(i,:));
    plot(X5(i,:),V(i,:));
end
hold off;
xlim([0 1]);
ylim([0 1]);
%}
%{
[X6, Y6] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
%a = 2 ./ (1 + Y4)
%b = Y4 .* a
%((b - a) / 2) * X.^2 + a .* X;
K = (((Y6 .* (2 ./ (1 + Y6))) - (2 ./ (1 + Y6))) ./ 2);% .* X6.^2;
L = (2 ./ (1 + Y6));% .* X6;
figure;
hold on;
for i = 1:size(Y6,1)
    plot(X6(i,:),L(i,:)./K(i,:));
end
T = X6.^2;
plot(X6(1,:),T,'g--o');
hold off;
xlim([0 1]);
ylim([0 1]);
% turns out the relation between cofA and cofB -> (z-1)/2 = cofA/cofB
% (cofB/cofA = 2/(z-1), and z = b/a (ratio of b to a)
% thus cofA = cofB * ((z-1) / 2)
% note that this is with pre-scaled a and b
%}

%{ 
% this seems to work for the CDF
[X7, Y7] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
S = ((Y7-1)./(Y7+1)) .* X7.^2;
D = (2./(Y7+1)) .* X7;
figure;
hold on;
for i = 1:size(Y7,1)
    plot(X7(i,:),S(i,:) + D(i,:));
end
hold off;
xlim([0 1]);
ylim([0 1]);
%}

%{
% and this for inverse CDF
[X7, Y7] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
S2 = (2./(Y7+1)) .* sqrt(X7);
D2 = ((Y7-1)./(Y7+1)) .* X7;
figure;
hold on;
for i = 1:size(Y7,1)
    plot(X7(i,:),S2(i,:) + D2(i,:));
end
hold off;
xlim([0 1]);
ylim([0 1]);
%}

%{
% combined, appears to not be an exact match to pure function calculation
%[X7, Y7] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
[X7, Y7] = meshgrid(0:0.005:1,[3]);
nth = 20;
S = ((Y7-1)./(Y7+1)) .* (X7.^nth); % these S formulas are equivelant
%S = (1-(2./(Y7+1))).* (X7.^nth);
D = (2./(Y7+1)) .* X7;
S2 = ((Y7-1)./(Y7+1)) .* ((X7).^(1/nth));
D2 = (2./(Y7+1)) .* X7;
%D2 = (1-(2./(Y7+1))) .* X7;
%S2 = (1-((2-Y7)./Y7)) .* ((X7).^(1/nth)) %trying inverse of ration
%D2 = ((2-Y7)./Y7) .* X7 % trying inverse of ratio

figure;
hold on;
for i = 1:size(Y7,1)
    plot(X7(i,:),S(i,:) + D(i,:));
    plot(X7(i,:),S2(i,:) + D2(i,:));
end
T = X.^nth;
R = X.^(1/nth);
plot(X(1,:),T,'g--o');
plot(X(1,:),R,'b--o');
hold off;
xlim([0 1]);
ylim([0 1]);
%}

%{
% WIP figuring out ratio inverse CDF
[X8, Y8] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
Z2 = (-1 + sqrt(1^2 - 4.*((Y8-1) ./ 2) .* (0 - (X8 .* 0.5.*(1+Y8))))) ./ (2 .* ((Y8-1) ./ 2));
%Q2 = ((((Y8-1) ./ 2).*(X8.^2) + 1.*X8 + 0)) ./ (0.5.*(1+Y8));
U = X8./Z2;
I = (X8.^2)./Z2;
figure;
hold on;
for i = 1:size(Y8,1)
    %plot(X8(i,:),U(i,:));
    %plot(X8(i,:),I(i,:));
    plot(X8(i,:),I(i,:)./U(i,:));
end
hold off;
xlim([0 1]);
ylim([0 1]);
%}

%{
% worked out n-formula for ratios scaled CDF
% combined, appears to not be an exact match to pure function calculation
%[X9, Y9] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
[X9, Y9] = meshgrid(0:0.005:1,[3]);
nth = 3;
S = ((Y9-1)./(nth+Y9-1)) .* (X9.^nth);
D = (nth ./ (nth+Y9-1)) .* X9;
SD = (((Y9-1)./nth) .* (X9.^nth) + X9) ./ (((Y9-1)./nth) + 1);
figure;
hold on;
for i = 1:size(Y9,1)
    plot(X9(i,:),S(i,:) + D(i,:));
    %plot(X9(i,:),SD(i,:));
end
T = X9.^nth;
R = X9.^(1/nth);
plot(X9(1,:),T,'g--o');
plot(X9(1,:),R,'b--o');
hold off;
xlim([0 1]);
ylim([0 1]);
legend
%}

[X9, Y9] = meshgrid(0:0.005:1,[1:0.2:2,3:0.5:5,6:1:10,20:10:100,200:100:1000,2000:1000:10000]);
%[X9, Y9] = meshgrid(0:0.005:1,[1 2 5 10 100]);
nth = 20;
S = ((Y9-1)./(nth+Y9-1)) .* (X9.^nth);
D = (nth ./ (nth+Y9-1)) .* X9;
Dtest = 1- (S+D);
SD = (((Y9-1)./nth) .* (X9.^nth) + X9) ./ (((Y9-1)./nth) + 1);
%S2 = ((Y9-1)./(nth+Y9-1)) .* (X9.^(1/nth));
%D2 = (nth ./ (nth+Y9-1)) .* X9;
%Y9 = 1 ./ Y9;
%S2 = ((-nth*Y9+Y9-1)./(Y9-1)) .* (X9.^(1/nth));
D2 = (1 ./ (nth ./ (nth+Y9-1))) .* X9;
S3 = ((1 ./ ((Y9-1)./(nth+Y9-1))) .* X9).^(1/nth);
%SD2 = ((((Y9.^nth).*X9-X9+1).^(1/nth)) - 1) ./ (Y9 -1);
SD3 = (-1 + (1^2 - 4.*((Y9-1) ./ nth) .* (0 - (X9 .* (1./nth).*((nth-1)+Y9)))).^(1./nth)) ./ (2 .* ((Y9-1) ./ nth));
Z = (-1 + sqrt(1^2 - 4.*((Y9-1) ./ 2) .* (0 - (X9 .* 0.5.*(1+Y9))))) ./ (2 .* ((Y9-1) ./ 2));
D4 = ((X9.*Y9).^nth);
figure;
hold on;
for i = 1:size(Y9,1)
    [col] = Rainbow(rand);
    %col = [min(max(0,cr),1) min(max(0,cr-1),1) min(max(0,cr-2),1)];
    %plot(X9(i,:),S(i,:) + D(i,:));
    %plot(X9(i,:),S(i,:));
    %plot(X9(i,:),D(i,:));
    %plot(X9(i,:),SD(i,:),'Color',col);
    %plot(X9(i,:),S2(i,:) + D2(i,:));
    %plot(X9(i,:),D2(i,:),'Color',col);
    %plot(X9(i,:),S3(i,:));
    %plot(X9(i,:),SD2(i,:),'Color',col);
    %plot(X9(i,:),SD3(i,:),'Color',col);
    %plot(X9(i,:),Z(i,:));
    %plot(X9(i,:),D4(i,:));
    plot(X9(i,:),Dtest(i,:));
end
T = X9.^nth;
R = X9.^(1/nth);
%plot(X9(1,:),T,'g--o');
%plot(X9(1,:),R,'b--o');
hold off;
xlim([0 1]);
ylim([0 1]);




function [Wi, bc, pdf] = WeightedTriangle(VertexValues)

% VertexValues = [A, B, C]
%
%       C
%      /\_ _ D
%     /  \
%   A/____\B

a = VertexValues(1);
b = VertexValues(2);
c = VertexValues(3);

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
%pU = 2 * acos(1 - rand) / pi;
pU = sqrt(rand);
%pU = rand;

if a == d
    ad = pU;
    %ad = sqrt(pU);
else
    ad = (-cofB2 + sqrt(cofB2^2 - 4*cofA2*(cofC2 - (pU*A)))) / (2*cofA2);
    %ad = ad^(1/1.7);
    %ad = sqrt(pU) * ad;
    %ad = sqrt(pU);
    %ad = sqrt(ad);
end


Wi = ad * adVector;

pdf = (dMinA*ad + a);

%th = -1/12*2*pi; %30 degrees clockwise
th = 2/12*2*pi; %60 degrees anti-clockwise
rm = [[cos(th) -sin(th)];[sin(th) cos(th)]];

Wi = transpose(rm*transpose(Wi));

end