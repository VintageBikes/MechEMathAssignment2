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