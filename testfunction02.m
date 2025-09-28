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



function [f_out,dfdx] = test_function02(X)
x1 = X(1);
x2 = X(2);
x3 = X(3);
y1 = 3*x1^2 + .5*x2^2 + 7*x3^2-100;
y2 = 9*x1-2*x2+6*x3;
f_out = [y1;y2];
dfdx = [6*x1,1*x2,14*x3;9,-2,6];
end