function [x, x_old, k] = n2_simpleiter(x0, H, g, epsilon)

x_old = x0;
x = H*x_old + g;
k = 1;

while abs(norm(x - x_old)) >= epsilon
    k = k + 1;
    x_old = x;

    x = H*x_old + g;
end

end
