function multinewtonsolver()
    X0 = rand(3,1);

    [~, J_analytical] = test_function02(X0);
    J_numerical = approximate_Jacobian01(@test_function02,X0);

    solver_params = struct();
    solver_params.dxmin = 1e-10;
    solver_params.ftol = 1e-10;
    solver_params.dxmax = 1e8;
    solver_params.max_iter = 200;
    solver_params.approx = 1;

    Xguess = randn(3,1);

    X_root = multivariate_Newton(@test_function02,Xguess,solver_params);

    disp(X_root);

    f_root = test_function02(X_root);

    disp(f_root);

    % disp(J_analytical)
    % disp(J_numerical)

end

function X_root = multivariate_Newton(fun,X,solver_params)
    dxmin = solver_params.dxmin;
    ftol = solver_params.ftol;
    dxmax = solver_params.dxmax;
    max_iter = solver_params.max_iter;
    approx = solver_params.approx;

    if approx
        fval = fun(X);
        J = approximate_Jacobian01(fun,X);
    else
        [fval,J] = fun(x);
    end

    delta_x = 1;

    count = 0;

    while count < max_iter && norm(delta_x) > dxmin && norm(fval) > ftol && norm(delta_x) < dxmax
        count = count+1;
        if approx
        fval = fun(X);
        J = approximate_Jacobian01(fun,X);
        else
        [fval,J] = fun(X);
        end
        
        delta_x = -J\fval;

        X = X + delta_x;


    end

    X_root = X;

end
function J = approximate_Jacobian01(fun,X)
    f0 = fun(X);

    J = zeros(length(f0),length(X));

    e_n = zeros(length(X),1);

    delta_X = 1e-6;

    for n = 1:length(X)
        e_n(n) = 1;

        f_plus = fun(X+e_n*delta_X);
        f_minus = fun(X-e_n*delta_X);

        J(:,n) = (f_plus-f_minus)/(2*delta_X);


        e_n(n) = 0;
    end 
end



function [f_out,dfdx] = test_function02(X)
x1 = X(1);
x2 = X(2);
x3 = X(3);
y1 = 3*x1^2 + .5*x2^2 + 7*x3^2-100;
y2 = 9*x1-2*x2+6*x3;
f_out = [y1;y2];
dfdx = [6*x1,1*x2,14*x3;9,-2,6];
end