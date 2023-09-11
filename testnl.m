% calculating trajectory of soccer ball kicked
% 4x4 nonlinear system
% z1 = x, z2 = y, z3 = theta, z4 = v

m = 0.45; g = 9.81;
rho = 1.29; s = 0.038;
W = 0;

theta = [pi/4 + 0.125 pi/4 pi/4 - 0.125 pi/4 - 0.25];
c = 0.2;

figure
fprintf('N   angle x   y  theta    speed    max(t)    min(t)   average(t)    W,   tn\n');
for a = 1:4
    z10 = 0; z20 = 0; z30 = theta(a); z40 = 22.5;
    z0 = [z10 z20 z30 z40]';
    dt = 0.001;
    maxit = 10; tol = 1e-12; % Newton's parameters

    z = z0; zi(a,1,:) = z; nit = [];
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
            if (norm(f, inf) <= tol), k = k-1; break, end
            J = [1 0 dt*(z(4)*sin(z(3))) -dt*cos(z(3));
                0 1 -dt*(z(4)*cos(z(3))) -dt*sin(z(3));
                0 0 1 - dt*g*sin(z(3))/z(4) -dt*g*cos(z(3))/(z(4)*z(4));
                0 0 dt*g*cos(z(3)) 1+dt*c*rho*s*z(4)/m];
            z = z - J\f(1:4,1);
        end
        zi(a, i+1,:) = z;
        nit(i) = k;
        if z(2) < 0
            x = zi(a, i-2:i,2);
            y = zi(a, i-2:i,1);
            h = [(i - 2)*dt, (i-1)*dt, i*dt];
            zi(a, i+1, 1) = interp1(x, y, 0, "cubic");
            tn = interp1(x, h, 0, "cubic");
            zi(a, i+1, 2) = 0;
            break
        end
        i = i + 1;
    end
    N1(a) = i + 1;
    fprintf('%4d %6.4f %8.5f %9.6f %8.5f %7.4f %3d %3d %5.2f %3.0f %7.4f\n', N1(a)-1, ...
        zi(a, 1, 3), zi(a, N1(a), 1), zi(a, N1(a), 2), zi(a, N1(a), 3), ...
        zi(a, N1(a), 4), max(nit), min(nit), sum(nit)/N1(a), W, tn);
end

%plot(zi(1, 1:200:N1(1)+1, 1), zi(1, 1:200:N1(1)+1, 2), 'r.', ...
     %zi(2, 1:200:N1(2)+1, 1), zi(2, 1:200:N1(2)+1, 2), 'b.', ...
     %zi(3, 1:200:N1(3)+1, 1), zi(3, 1:200:N1(3)+1, 2), 'g.', ...
     %zi(4, 1:200:N1(4)+1, 1), zi(4, 1:200:N1(4)+1, 2), 'k.');
%legend('\theta = \pi/4 + 0.125', ...
%  '\theta = \pi/4', ...
%  '\theta = \pi/4 - 0.125', ...
%  '\theta = \pi/4 - 0.25');
%axis equal
%title('simulating movement of ball');
%xlabel('horizontal position');
%ylabel('vertical position');
