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
