clear ; close all; clc

A = [3.278164  1.046583 -1.378574; 
     1.046583  2.975937  0.934251; 
    -1.378574  0.934251  4.836173]
%A = [-0.81417 -0.01937 0.41372;
%     -0.01937 0.54414 0.00590;
%     0.41372 0.00590 -0.81445]
%eig(A)

epsilon = 1e-6

[lmbd, V] = eigenjacob(length(A), A, epsilon)

disp('')
disp('================================')
disp('')
disp('Решение методом Якоби: ')
disp('Собственные числа')
disp(lmbd)
disp('Собственные векторы')
disp(V)
disp('')
disp('================================')
%[u,v] = eig(A)

epsilon = 0.001;
[lmbd, v, k] = eigenpow(length(A), A, epsilon);

disp('')
fprintf('Степенной метод, наибольшее собственное число: %d\n', lmbd)
fprintf('Собственный вектор (его норма равна %d):\n', norm(v, 2))
disp(v)
fprintf('Количество шагов: %d\n', k)
fprintf('Невязка: %d\n', norm(A*v - lmbd*v))
disp('')
disp('================================')

epsilon = 0.000001;
[lmbd, v, k] = eigenscal(length(A), A, epsilon);

disp('')
fprintf('Метод скалярных произведений, наибольшее собственное число: %d\n', lmbd)
fprintf('Собственный вектор (его норма равна %d):\n', norm(v, 2))
disp(v)
fprintf('Количество шагов: %d\n', k)
fprintf('Невязка: %d\n', norm(A*v - lmbd*v))
disp('')
disp('================================')

epsilon = 0.001;
lmbd_max = lmbd;
v_max = v;
B = A - lmbd_max * eye(length(A));

[lmbd, v, k] = eigenpow(length(B), B, epsilon);
lmbd = lmbd + lmbd_max;
%v = v + v_max;

disp('')
fprintf('Противоположная граница спектра (степенным методом)\n')
fprintf('Собственное число: %d\n', lmbd)
fprintf('Собственный вектор (его норма равна %d):\n', norm(v, 2))
disp(v)
fprintf('Количество шагов: %d\n', k)
fprintf('Невязка: %d\n', norm(A*v - lmbd*v))
disp('')
disp('================================')

epsilon = 0.001;
[lmbd, v, k] = wielandt(length(A), A, epsilon, lmbd_max);


disp('')
fprintf('Метод Виландта, наибольшее по модулю собственное число:\n')
fprintf('Собственное число: %.15d\n', lmbd)
fprintf('Собственный вектор (его норма равна %d):\n', norm(v, 2))
disp(v)
fprintf('Количество шагов: %d\n', k)
fprintf('Невязка: %d\n', norm(A*v - lmbd*v))
disp('')

