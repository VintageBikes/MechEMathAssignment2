function V_p = projectile_traj(X)
    theta = X(1);
    t = X(2);
    g = 2.3; %gravity in m/sˆ2
    v0 = 14; %initial speed in m/s
    px0 = 2; %initial x position
    py0 = 4; %initial y position
    
    %compute position vector
    V_p = [v0*cos(theta)*t+px0; -.5*g*t.^2+v0*sin(theta)*t+py0];
end