X = [1; 2; 3];
[F, J] = test_function01(X)

%the function name and input/output variable names
%are just what I chose, you can use whatever names you'd like
function [F, J] = test_function01(X)
    f1 = X(1)^2 + X(2)^2 - 6 - X(3)^5;
    f2 = X(1)*X(3) + X(2) - 12;
    f3 = sin(X(1) + X(2) + X(3));
    F = [f1; f2; f3];

    pd_x1 = [
        2*X(1) + X(2)^2 - 6 - X(3)^5
        X(3) + X(2) - 12
        cos(X(1) + X(2) + X(3))
    ];
    pd_x2 = [
        X(1)^2 + 2*X(2) - 6 - X(3)^5
        X(1)*X(3) - 11
        cos(X(1) + X(2) + X(3))
    ];
    pd_x3 = [
        X(1)^2 + X(2)^2 - 6 - 5*X(3)^4
        X(1) + X(2) - 12
        cos(X(1) + X(2) + X(3))
    ];
    J = [pd_x1, pd_x2, pd_x3];
end