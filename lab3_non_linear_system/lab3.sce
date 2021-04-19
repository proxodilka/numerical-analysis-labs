
function [res]=f1(x)
    res = sin(x(1, :) + x(2, :)) - 1.2*x(1, :) - 0.1
endfunction

function [res]=f2(x)
    res = x(1, :)^2 + x(2, :)^2 - 1
endfunction

function [res]=f(x)
    res = [f1(x); f2(x)]
endfunction

function [res]=jacob_exact(x)
    // "Аналитический" Якобиан функции
    res = [
        cos(x(1, :) + x(2, :)) - 1.2,    cos(x(1, :))*cos(x(2, :)) - sin(x(1, :))*sin(x(2, :));
        2 * x(1, :)                 ,    2 * x(2, :)
    ]
endfunction

function [res]=jacob_approx(x)
    // Якобиан функции, посчитанный численными методами
    res = numderivative(f, x)
endfunction

function [x0]=newtoon_method(f, x0, eps, jacob, max_iter, verbose)
    if ~exists("eps") then eps=10e-6 end
    if ~exists("max_iter") then max_iter=100 end
    if ~exists("verbose") then verbose = %F end
    
    for i=1:max_iter
        x_prev = x0
        x0 = x0 - jacob(x0)^-1 * f(x0)
        
        err = norm(x0 - x_prev)
        if verbose then
            printf("\nШаг %i:\n--------", i)
            disp(["x:", "f(x):"; string(x0), string(f(x0))])
        end
        if err < eps then break end
    end
endfunction

function [res1, res2, res3]=verbose_solver(x0, eps)
    disp("f(x) в начальном приближении:", f(x0))
    disp("Якобиан в начальном приближении:", jacob_exact(x0))
    disp("Обратный Якобиан в начальном приближении:", jacob_exact(x0)^-1)
    printf("\n--------------- Метод Ньютона с ''приближенным'' Якобианом ---------------\n")
    res1 = newtoon_method(f, x0, eps=eps, jacob=jacob_approx, verbose=%T)
    printf("\n--------------- Метод Ньютона с точным Якобианом ---------------\n")
    res2 = newtoon_method(f, x0, eps=eps, jacob=jacob_exact, verbose=%T)
    printf("\n---------------------  SciLab.fsolve  ------------------------\n")
    res3 = fsolve(x0, f, fjac=jacob_exact, tol=eps)
    disp(["x:", "f(x):"; string(res3), string(f(res3))])
endfunction

// 
x = [-1:0.01:1]'
plot2d(x, [asin(1.2*x + 0.1) - x, sqrt(1 - x^2), -sqrt(1 - x^2)], [1, 2, 2])
legend("sin(x + y) - 1.2x - 0.1 = 0", "x^2 + y^2 - 1 = 0")

// Визуально можно определить два корня:
x01 = [-0.89, -0.45]'
x02 = [0.74, 0.67]'

printf("------------ Поиск первого решения системы ------------")
[x11, x12, x13] = verbose_solver(x01, eps=10e-6)
printf("------------ Поиск второго решения системы ------------")
[x21, x22, x23] = verbose_solver(x02, eps=10e-6)
