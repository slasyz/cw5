function [] = main2disp(y, y_exact)

h_old = 0.1;
h_new = 0.001;

for i = 1:5
    disp(abs(y_exact(i*100, :) - y(i*100, :)));
end

end
