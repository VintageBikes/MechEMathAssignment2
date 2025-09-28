function projectile_solver()
    solver_params = struct();
    solver_params.dxmin = 1e-10;
    solver_params.ftol = 1e-10;
    solver_params.dxmax = 1e8;
    solver_params.max_iter = 200;
    solver_params.approx = 1;

    Xguess = [pi/4;11] %Initial guess: 45 degrees and 2 seconds

    X_root = multivariate_Newton(@collision_function,Xguess,solver_params);

    disp(X_root);

    f_root = collision_function(X_root);

    disp(f_root);

    run_simulation(X_root(1), X_root(2));

end

