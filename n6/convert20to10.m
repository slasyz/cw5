function [y] = convert20to10(y_orig)

y = zeros(11, 1);
for i = 0:10
    %fprintf('%d <- %d\n', i, 2*i);
    y(i+1) = y_orig(2*i + 1);
end
