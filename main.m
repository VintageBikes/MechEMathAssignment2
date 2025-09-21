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

%% Part 3
J = approximate_jacobian(@test_function01, X)
test_numerical_jacobian()

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

%%
%Implementation of finite difference approximation
%for Jacobian of multidimensional function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%J: approximation of Jacobian of fun at x
function J = approximate_jacobian(fun, X)
    h = 1e-6;
    J = zeros(length(fun(X)), length(X));
    for j = 1:length(X)
        e = zeros(length(X), 1);
        e(j) = 1;
        J(:,j) = (fun(X + h*e) - fun(X - h*e)) / (2*h);
    end
end

%function to test numerical differentiation function
function test_numerical_jacobian()
    %number of tests to perform
    num_tests = 100;
    %iterate num_tests times
    for n = 1:num_tests
        %generate a randomized input and output dimension
        input_dim = randi([1,15]);
        output_dim = randi([1,15]);
        %generate a input_dim x input_dim x output_dim matrix stack A
        A = randn(input_dim,input_dim,output_dim);
        %generate a matrix, B of dimension output_dim x input_dim
        B = randn(output_dim,input_dim);
        %generate a column vector, C of height output_dim
        C = randn(output_dim,1);
        %create a new test function
        %this is essentially a random second-order (quadratic) function
        %with input dimension input_dim and output dimension output_dim
        test_fun = @(X) jacobian_test_function(X,A,B,C);
        X_guess = randn(input_dim,1);
        %evaluate numerical Jacobian of test_fun
        %use whatever your function name was here!
        J_numerical = approximate_jacobian(test_fun,X_guess);
        %compute the analytical jacobian of jacobian_test_function
        J_analytical = B;
        for n = 1:output_dim
            J_analytical(n,:)=J_analytical(n,:)+X_guess'*A(:,:,n);
            J_analytical(n,:)=J_analytical(n,:)+X_guess'*A(:,:,n)';
        end
        %compare with Jacobian of A
        largest_error = max(max(abs(J_numerical-J_analytical)));
        %if J is not close to A, print fail.
        if largest_error>1e-7
            disp('fail!');
        end
    end
end

%computes a quadratic function on input X
function f_val = jacobian_test_function(X,A,B,C)
    output_length = length(C);
    f_val = B*X+C;
    for n = 1:output_length
        f_val(n)=f_val(n)+(X'*A(:,:,n)*X);
    end
end