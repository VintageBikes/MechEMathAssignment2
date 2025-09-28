function F = collision_function(X)
    theta = X(1);
    t = X(2);
    F = projectile_traj([theta; t]) - target_traj(t);
end
