function [x,flag,iter] = secant(fun, x1, x2, tol, maxit)

if nargin < 5
    maxit = 100;
end

f1 = fun(x1);
if abs(f1) < tol
    x = x1;
    return;
end

f2 = fun(x2);
if abs(f2) < tol
    x = x2;
    return;
end

flag = 1;

for iter = 1:maxit
    if abs(f1-f2) < 1e-12
        flag = 0;
        break;
    end
    x = x1 - (x1-x2)/(f1-f2) * f1;
    f = fun(x);
    
    fprintf('%d\t%.3e\t%.3e\n', iter, x, abs(f));
    
    if abs(f) < tol
        flag = 0;
        break;
    else
        x2 = x1;
        f2 = f1;
        x1 = x;
        f1 = f;
    end
end

end