function [lmbd, X] = eigenjacob(n, A, epsilon)

X = eye(n);

while true
    [M, I] = max(abs(triu(A, 1)(:)));
    [i, j] = ind2sub(size(A), I);

    if M < epsilon
        break
    end

    d = sqrt((A(i,i) - A(j,j))^2 + 4*A(i,j)^2);
    c = sqrt(1/2 * (1 + abs(A(i,i) - A(j,j))/d));
    s = sign(A(i,j)*(A(i,i) - A(j,j))) * sqrt(1/2 * (1 - abs(A(i,i) - A(j,j))/d));

    %phi = atan((2 * A(i, j)) / (A(i, i) - A(j, j))) / 2;
    %c = cos(phi);
    %s = sin(phi);

    A_new = A;
    i_k = i; j_k = j;

    for i = 1:n
        A_new(i, i_k) = A_new(i_k, i) = c*A(i, i_k) + s*A(i, j_k);
        A_new(i, j_k) = A_new(j_k, i) = -s*A(i, i_k) + c*A(i, j_k);
    end
    A_new(i_k, i_k) = c^2 * A(i_k, i_k) + 2*c*s*A(i_k, j_k) + s^2*A(j_k, j_k);
    A_new(j_k, j_k) = s^2 * A(i_k, i_k) - 2*c*s*A(i_k, j_k) + c^2*A(j_k, j_k);
    A_new(i_k, j_k) = A_new(j_k, i_k) = (c^2 - s^2)*A(i_k, j_k) + c*s*(A(j_k, j_k) - A(i_k, i_k));

    A = A_new;

    for i = 1:n
        temp = X(i, i_k);
        X(i,i_k) = c*X(i, i_k) + s*X(i, j_k);
        X(i,j_k) = -s*temp + c*X(i, j_k);
    end;
end

lmbd = diag(A);

end
