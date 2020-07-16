
%{
    [I]=G________F  +y
	    /|     /|    ^
	   / |    / |    |
	 H/__|__E/  |
	  |  C--|--/B       +z
	  | /   | /        /
	  |/____|/       -z
	 D     A=(0,0,0)
	    +x<->-x
%}

%max value = 101
%min value = 0

% situation values where floodfill enterted from negative z direction
% and left side is completly opaque
datapoints = [
    0, ...      %A  green
    1, ...    %B    magenta
    101, ...    %C  cyan
    101, ...    %D  red
    0, ...    %E    orange
    1, ...    %F    blue
    101, ...    %G  pink
    101];       %H  black
%diagonal backward outerpoint
datapoints = [
    0, ...      %A  green
    1, ...    %B    magenta
    1, ...    %C  cyan
    101, ...    %D  red
    0, ...    %E    orange
    1, ...    %F    blue
    1, ...    %G  pink
    101];       %H  black
%diagonal forward outerpoint
datapoints = [
    0, ...      %A  green
    101, ...    %B    magenta
    1, ...    %C  cyan
    0, ...    %D  red
    0, ...    %E    orange
    101, ...    %F    blue
    1, ...    %G  pink
    0];       %H  black
%diagonal backward innerpoint
datapoints = [
    104, ...      %A  green
    4, ...    %B    magenta
    105, ...    %C  cyan
    105, ...    %D  red
    104, ...    %E    orange
    4, ...    %F    blue
    105, ...    %G  pink
    105];       %H  black
%diagonal forward innerpoint
datapoints = [
    106, ...      %A  green
    107, ...    %B    magenta
    107, ...    %C  cyan
    6, ...    %D  red
    106, ...    %E    orange
    107, ...    %F    blue
    107, ...    %G  pink
    6];       %H  black
n = 10000;
Kd = [1, 1, 1];

cumalative = [0, 0, 0, 0, 0, 0, 0, 0];
dirs = cell(8,1);
for i = 1:n
       [F, Wi, Pdf, Face] = OctoGradient(Kd, datapoints);
       for f = 1:8
           if Face(f) == 1
               dirs{f} = [dirs{f}; Wi];
               break;
           end
       end
       cumalative = cumalative + Face;
       if rem(i, n / 10) == 0
           i / n;
       end
end

figure;
colors = {'g','m','c','r',[1,0.5,0],'b',[1,0.5,0.75],'k'};
holding = 0;
for c = 1:8
    if size(dirs{c}, 1) > 0
        scatter3(dirs{c}(:,1), dirs{c}(:,2), dirs{c}(:,3), 5, colors{c});
        if not(holding)
            hold on
            holding = 1;
        end
    end
end
hold off
xlabel('X')
ylabel('Y')
zlabel('Z')

%plot total directions
figure;
plot(1:size(cumalative,2),cumalative);

%Plot chances
figure;
plot(1:size(cumalative,2),cumalative/n);