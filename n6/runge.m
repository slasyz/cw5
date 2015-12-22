function [y] = runge(y_orig, y_orig_small)

y = y_orig_small;

%for i = 0:10
%    y(i+1) = 2*y_orig(2*i+1) - y_orig(i+1);
%end

for i = 1:11
    %y(i) = 2 * y_orig(2*i - 1) - y_orig_small(i);
    y(i) = y_orig(2*i - 1) + (y_orig(2*i - 1) - y_orig_small(i))/(2^2-1);
end;
