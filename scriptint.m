a = 0; b = 1; M = 1000; v = linspace(a, b, M);
u1 = v.^(5/2); u2 = v.^(3/2); u3 = v.^(1/2);
fprintf('   n  v.^(5/2)  v.^(3/2)  v.^(1/2)  respective orders\n');
for nn = 1:20
    n(nn) = 2^(nn+4);
    x = linspace(a, b, n(nn)+1);
    y1 = x.^(5/2); y2 = x.^(3/2); y3 = x.^(1/2);
    e1(nn) = max(abs(interp1(x, y1, v, 'linear') - u1));
    e2(nn) = max(abs(interp1(x, y2, v, 'linear') - u2));
    e3(nn) = max(abs(interp1(x, y3, v, 'linear') - u3));
    fprintf('%4d %9.2e %9.2e %9.2e  ', n(nn), e1(nn), e2(nn), e3(nn));
    if nn > 1
         c1(nn) = log(e1(nn-1)/e1(nn))/log(2);
         c2(nn) = log(e2(nn-1)/e2(nn))/log(2);
         c3(nn) = log(e3(nn-1)/e3(nn))/log(2);
         fprintf('%5.2f %5.2f %5.2f', c1(nn), c2(nn), c3(nn));
    end
    fprintf('\n');
end
figure;
loglog(n, e1, 'r-', n, e2, 'g-', n, e3, 'b-', 'LineWidth', 2);
xlabel('number of multiplications');
ylabel('Error for different functions');
title('testing for convergence for linear spline interpolation');
legend('y1 = x^{5/2}', 'y2 = x^{3/2}', 'y3(x) = x^{1/2}');
grid on;
