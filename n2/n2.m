clear ; close all; clc

%A = input('Введите A: ')
A = [3.278164  1.046583 -1.378574; 
     1.046583  2.975937  0.934251; 
    -1.378574  0.934251  4.836173]
n = length(A)
%b = input('Введите b: ')
b = [-0.527466; 2.526877; 5.165441]
epsilon = 1e-3
x0 = zeros(n, 1)

x_gauss = A\b
disp('================================')
disp('Решение методом Гаусса: ')
disp(x_gauss)
disp('================================')
disp('')

%H = zeros(n)
%for i = 1:n
%    for j = 1:n
%        if i == j
%            H(i,j) = 0;
%        else
%            H(i,j) = - A(i,j) / A(i,i);
%        end
%    end
%    g(i) = b(i) / A(i,i);
%end

%disp('Старенькое: ')
%H
%g

% =======================
% Приводим к нужному виду
% =======================
D = zeros(n);
for i = 1:n
    D(i,i) = A(i,i);
end

H = eye(n) - inv(D)*A;
g = inv(D)*b;

H
g

fprintf('Норма матрицы H: %f\n', norm(H, Inf));

%k_aprior = (log(epsilon) - log(norm(x0)) + log(1-norm(H))) / log(norm(H)) - 1
k_aprior = log(epsilon*(norm(x0) + norm(g)/(1-norm(H)))^(-1))/log(norm(H));
fprintf('Априорная оценка k: %f\n', k_aprior)

% ======================
% Метод простой итерации
% ======================

[x, x_old, k] = n2_simpleiter(x0, H, g, epsilon);
x_lust = x_old + (x-x_old)/(1-max(abs(eig(H))));

disp('==============================')
disp('Метод простой итерации, ответ: ')
disp(x)
fprintf('Найдено за %d шагов\n', k)
fprintf('Фактическая погрешность: %d\n', norm(x-x_gauss))
fprintf('Апостериорная оценка: %f\n', norm(H)*norm(x-x_old)/(1-norm(H)))
fprintf('Априорная оценка: %f\n', norm(H)*norm(x0) + norm(H)^k*norm(g)/(1-norm(H)))
disp('')
disp('Уточнение по Люстернику: ')
disp(x_lust)
fprintf('Фактическая погрешность: %d\n', norm(x_lust-x_gauss))
disp('==============================')
disp('')

% =============
% Метод Зейделя
% =============

H_l = tril(H); H_r = triu(H);
H_seid = (eye(n) - H_l)^(-1) * H_r;

[x, x_old, k] = n2_seid(x0, H, g, epsilon);

disp('=================================')
disp('Метод Зейделя, ответ: ')
disp(x)
fprintf('Найдено за %d шагов\n', k)
fprintf('Фактическая погрешность: %d\n', norm(x-x_gauss))
fprintf('Апостериорная оценка: %f\n', norm(H)*norm(x-x_old)/(1-norm(H)))
%fprintf('Априорная оценка: %f\n', norm(H)*norm(x0) + norm(H)^k*norm(g)/(1-norm(H)))
disp('')
fprintf('Спектральный радиус H: %f\n', max(abs(eig(H))))
fprintf('Спектральный радиус H_seid: %f\n', max(abs(eig(H_seid))))
disp('=================================')
disp('')

% ========================
% Метод верхней релаксации
% ========================

po = max(abs(eig(H)))
q = 2 / (1 + sqrt(1 - po^2))
%q = 1
[x_opt,   x_old, k] = n2_upprel(x0, H, g, epsilon, q);
[x_minus, x_old, k] = n2_upprel(x0, H, g, epsilon, q-0.2);
[x_plus,  x_old, k] = n2_upprel(x0, H, g, epsilon, q+0.2);

disp('================================')
disp('Метод верхней релаксации, ответ:')
fprintf('Найдено за %d шагов\n', k)
fprintf('Фактическая погрешность при оптимальном q: %.16f\n', norm(x_opt-x_gauss))
fprintf('Фактическая погрешность при q-0.1: %.16f\n', norm(x_minus-x_gauss))
fprintf('Фактическая погрешность при q+0.1: %.16f\n', norm(x_plus-x_gauss))
disp('================================')
disp('')
