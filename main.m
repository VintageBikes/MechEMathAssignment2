X = [1; 2; 3];

%% Part 1
[F, J] = test_function01(X)

%% Part 2
dxtol = 1e-14;
ytol = 1e-14;
max_iter = 200;
dfdxmin = 1e-8;

x_root = multidimensional_newton_method(@test_function01, X, dxtol, ytol, max_iter, dfdxmin)
[F, J] = test_function01(x_root)

%%
function [F, J] = test_function01(X)
    f1 = X(1)^2 + X(2)^2 - 6 - X(3)^5;
    f2 = X(1)*X(3) + X(2) - 12;
    f3 = sin(X(1) + X(2) + X(3));
    F = [f1; f2; f3];

    pd_x1 = [
        2*X(1)
        X(3)
        cos(X(1) + X(2) + X(3))
    ];
    pd_x2 = [
        2*X(2)
        1
        cos(X(1) + X(2) + X(3))
    ];
    pd_x3 = [
        -5*X(3)^4
        X(1)
        cos(X(1) + X(2) + X(3))
    ];
    J = [pd_x1, pd_x2, pd_x3];
end

%%
function x_root = multidimensional_newton_method(fun, X_guess, dxtol, ytol, max_iter, dfdxmin)
    delta_x = 2 * dxtol;
    [F_X, J_X] = fun(X_guess);

    count = 0;
    while count < max_iter && max(abs(delta_x)) > dxtol && max(abs(F_X)) > ytol && max(abs(J_X(:))) > dfdxmin
        count = count + 1;
        delta_x = J_X \ -F_X;
        X_guess = X_guess + delta_x;
        [F_X, J_X] = fun(X_guess);
    end
    x_root = X_guess;
end