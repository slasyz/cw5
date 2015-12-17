function [y] = adams(y_exact, h, A)

n = length(A);
N = 0.5/h; % количество точек в сетке
y = zeros(N, 2);

W = inv(eye(n) - 5*h*A/12);
W1 = W * (eye(n) + 2*h*A/3);
W2 = W * (h*A/12);

% Находим решения через два предыдущих
y(1, :) = y_exact(1, :);
y(2, :) = y_exact(2, :);

for i = 3:N
    y(i, :) = W1 * y(i - 1, :)' + W2 * y(i - 2, :)';
end

end
