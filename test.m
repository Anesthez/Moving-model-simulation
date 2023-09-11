% calculating trajectory of soccer ball kicked
% 4x4 nonlinear system
% z1 = x, z2 = y, z3 = theta, z4 = v

m = 0.45; g = 9.81;
rho = 1.29; s = 0.038;
c = 0; W = 0; % c = 0.2;

theta = pi/4 - 0.25;

fprintf('it      x                y             theta            v             resid\n');
%c = 0.2;
z10 = 0; z20 = 0; z30 = theta; z40 = 22.5;
z0 = [z10 z20 z30 z40]';
dt = 0.001;
maxit = 10; tol = 1e-12; % Newton's parameters

z = z0; zi(:, 1) = z; nit = [];
i = 1;
while(true) % for (c)-(d)-(e) change this
    zinit = z; % initial guess for Newton's of the i-th step
    tn = i*dt; % new timestep
    W = 0; % may be changed depending on time and subquestion
    for k = 1:maxit
        f(1, 1) = z(1) - zinit(1) - dt*(z(4)*cos(z(3)) + W);
        f(2, 1) = z(2) - zinit(2) - dt*(z(4)*sin(z(3)));
        f(3, 1) = z(3) - zinit(3) + dt*g*cos(z(3))/z(4);
        f(4, 1) = z(4) - zinit(4) + dt*(c*rho*s*z(4)*z(4)/(2*m) + g*sin(z(3)));
        fprintf('%2d %16.14f %16.14f %16.14f %15.12f %8.2e\n', k-1, z(1:4),norm(f,inf));
        if (norm(f, inf) <= tol), k = k-1; break, end
        J = [1 0 dt*(z(4)*sin(z(3))) -dt*cos(z(3));
                0 1 -dt*(z(4)*cos(z(3))) -dt*sin(z(3));
                0 0 1 - dt*g*sin(z(3))/z(4) -dt*g*cos(z(3))/(z(4)*z(4));
                0 0 dt*g*cos(z(3)) 1+dt*c*rho*s*z(4)/m];
        z = z - J\f(1:4,1);
    end
    zi(:, i+1) = z;
    nit(i) = k;
    if z(2) < 0, break, end
    i = i + 1;
end
N1 = i+1;
fprintf('%4d %6.4f %8.5f %9.6f %8.5f %7.4f %3d %3d %5.2f %3.0f\n', N1-1, ...
        zi(3, 1), zi(1, N1), zi(2, N1), zi(3, N1), zi(4, N1), ...
        max(nit), min(nit), sum(nit)/N1, W);

