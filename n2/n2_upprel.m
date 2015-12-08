function [x, x_old, k] = n2_upprel(x0, H, g, epsilon, q)

x_old = x0;
x = n2_upprel_step(x_old, H, g, q);
k = 1;

while abs(norm(x - x_old)) >= epsilon
    k = k + 1;
    x_old = x;
    x = n2_upprel_step(x_old, H, g, q);
end

end

function [x] = n2_upprel_step(x_old, H, g, q)
    n = length(x_old)
    x = x_old

    for i = 1:n
        s = 0;
        for j = 1:i-1
            s = s + H(i,j)*x(j);
        end
        for j = i+1:n
            s = s + H(i,j)*x_old(j);
        end
        s = s - x_old(i) + g(i);

        x(i) = x_old(i) + q*s;
    end
end
