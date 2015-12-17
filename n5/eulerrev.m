function [y] = eulerrev(y0, h, A)

n = length(A);
N = 0.5/h; % количество точек в сетке
y = zeros(N, 2);

W = inv(eye(n) - h*A);

% Находим решения через одно предыдущее
y(1, :) = W * y0;
for i = 2:N
    y(i, :) = W * y(i - 1, :)';
end

end
