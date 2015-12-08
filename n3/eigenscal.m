function [lmbd, v, k] = eigenscal(n, A, epsilon)

k = 0;
%Y0 = rand(n,1);
Y0 = eye(n, 1);

while true
    k = k + 1;

    Y = A*Y0;
    if k > 1
        lmbd_old = lmbd;
    end
    lmbd = (Y' * Y0) / (Y0' * Y0);
    Y0 = Y / norm(Y);

    if k > 1
        if abs(lmbd - lmbd_old) <= epsilon
            break
        end 
    end 
end

%[v, D] = eig(A);
%v = v(:, 1);
v = Y0;

end
