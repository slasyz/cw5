function [lmbd, v, k] = wielandt(n, A, epsilon, lmbd0)

k = 0;
%Y0 = rand(n,1);
Y0 = eye(n, 1);
lmbd = lmbd0;

while true
    k = k + 1;

    if k > 1
        lmbd_old = lmbd;
    end 

    W = A - lmbd * eye(n);
    %lmbd
    %W
    Y = W \ Y0;
    Y0 = Y / norm(Y, 2);

    [mu, v, _] = eigenscal(n, inv(W), epsilon);
    %mu
    lmbd = 1/mu + lmbd;

    if k > 1
        if abs(lmbd - lmbd_old) <= epsilon
            break
        end 
    end 
end

%[v, D] = eig(A);
%v = v(:, 1);
%v = Y0;

end
