function [y] = approx(alpha1, beta1, alpha2, beta2, alpha, beta, a, b, n)

p = @(x) -(4-x)/(5-2*x);
q = @(x) (1-x)/2;
r = @(x) log(3+x)/2;
f = @(x) 1 + x/3;

h = (b-a)/n;

x = A = B = C = G = zeros(n+1, 1);
y = zeros(n+1, 1);
for i = 0:n
    x(i+1) = a + h*i;
end;

for i=2:n
    A(i) = p(x(i))/h^2 - q(x(i))/(2*h);
    B(i) = -r(x(i)) + 2*p(x(i))/h^2;
    C(i) = p(x(i))/h^2 + q(x(i))/(2*h);
    G(i) = f(x(i));
end

G(1) = alpha;
A(1) = 0;
B(1) = -(alpha1+alpha2/h);
C(1) = -alpha2/h;
G(n+1) = beta;
B(n+1) = -(beta1+beta2/h);
A(n+1) = -beta2/h;
C(n+1) = 0;

s = zeros(n+1, 1);
t = zeros(n+1, 1);

s(1) = C(1) / B(1);
t(1) = - G(1) / B(1);

for i = 2:n+1
    s(i) = C(i)/(B(i) - A(i)*s(i-1));
    t(i) = (A(i)*t(i-1) - G(i)) / (B(i) - A(i)*s(i-1));
end;

y(n+1) = t(n+1);
for i = n:-1:1
    y(i)=y(i+1)*s(i)+t(i);
end;
