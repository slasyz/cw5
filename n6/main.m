clear ; close all; clc;

format long;

p = @(x) -(4-x)/(5-2*x);
q = @(x) (1-x)/2;
r = @(x) log(3+x)/2;
f = @(x) 1 + x/3;

a = -1;
b = 1;
alpha1 = beta1 = 1;
alpha2 = beta2 = 0;
alpha = beta = 0;

x = zeros(11, 1);
for i = 0:10
    x(i+1) = -1 + 0.2 * i;
end

y_exact = [0.0000000000; 0.1175712552; 0.2254837911; 0.3173469344; 0.3865809340;
           0.4266772775; 0.4314721840; 0.3954079821; 0.3137602247; 0.1828160113; 0.0000000000];

% Мда
y_approx_10 = approx(alpha1, beta1, alpha2, beta2, alpha, beta, a, b, 10);
y_approx_20 = approx(alpha1, beta1, alpha2, beta2, alpha, beta, a, b, 20);
y_approx_20_small = convert20to10(y_approx_20);
y_runge = runge(y_approx_20, y_approx_20_small);

r1 = abs(y_exact - y_approx_10);
r2 = abs(y_exact - y_approx_20_small);
r3 = abs(y_exact - y_runge);

A_head = '      x        Y_ex     Y n=10   |Y-Yex|    Y n=20   |Y-Yex|     Y_r     |Y-Yex|';
A = [x y_exact y_approx_10 r1 y_approx_20_small r2 y_runge r3];

disp(A_head)
disp(A)
