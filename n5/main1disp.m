function [] = main1disp(y, y_exact)

for i = 1:5
    disp(abs(y_exact(i, :) - y(i, :)));
end

end
