datapoints = [0,101,101,101,101,101,101,101];
n = 100000000;
Kd = [1, 1, 1];

cumalative = [0, 0, 0, 0, 0, 0, 0, 0];
%dirs = zeros(n,3);
for i = 1:n
       [F, Wi, Pdf, Face] = OctoGradient(Kd, datapoints);
       %dirs(i,:) = Wi;
       cumalative = cumalative + Face;
       if rem(i, 10000000) == 0
           i / n
       end
end

%figure;
%scatter3(dirs(:,1), dirs(:,2), dirs(:,3));
%figure;
%plot(1:size(cumalative,2),cumalative);
figure;
plot(1:size(cumalative,2),cumalative/n);