%Implementation of finite difference approximation
%for Jacobian of multidimensional function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%J: approximation of Jacobian of fun at x
function J = approximate_jacobian(fun, x)
    f0 = fun(x);
    J = zeros(length(f0),length(x));
    e_n = zeros(length(x),1);
    delta_X = 1e-6;

    for n = 1:length(x)
        e_n(n) = 1;
        f_plus = fun(x+e_n*delta_X);
        f_minus = fun(x-e_n*delta_X);
        J(:,n) = (f_plus-f_minus)/(2*delta_X);
        e_n(n) = 0;
    end 
end