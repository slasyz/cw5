function [x, x_old, k] = n2_seid(x0, H, g, epsilon)

x_old = x0;
x = n2_seid_step(x_old, H, g);
k = 1;

while abs(norm(x - x_old)) >= epsilon
    k = k + 1;
    x_old = x;
    x = n2_seid_step(x_old, H, g);
end

end

function [x] = n2_seid_step(x_old, H, g)
    n = length(x_old);
    x = x_old;

    for i = 1:n
        s = 0;
        for j = 1:i-1
            s = s + H(i,j)*x(j);
        end
        for j = i:n
            s = s + H(i,j)*x_old(j);
        end
        s = s + g(i);

        x(i) = s;
    end

end
