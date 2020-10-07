function [Color] = Rainbow(r)
cr = r * 4;
if cr <= 1
    col = [1 cr 0];
elseif cr <= 2
    cr = cr - 1;
    col = [(1-cr) 1 0];
elseif cr <= 3
    cr = cr - 2;
    col = [0 1 cr];
else %cr <= 4
    cr = cr - 3;
    col = [0 (1-cr) 1];
end
Color = col;
end