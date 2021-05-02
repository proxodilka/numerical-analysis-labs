
function [res]=f1(x)
    res = sin(x(1, :) + x(2, :)) - x(2, :) - 1.5
endfunction

function [res]=f2(x)
    res = x(1, :) + cos(x(2, :) - 0.5) - 0.5
endfunction

function [res]=f(x)
    res = [f1(x); f2(x)]
endfunction

function [res]=jacob_exact(x)
    // "Аналитический" Якобиан функции
    res = [
        cos(x(1, :) + x(2, :)),    cos(x(1, :) + x(2, :)) - 1;
        1                     ,    sin(0.5 - x(2, :))
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

function [res]=f1_p1(y, n)
    n = ones(y) * n
    res = 2*%pi*n - y + asin(y + 1.5)
endfunction

function [res]=f1_p2(y, n)
    n = ones(y) * n
    res = 2*%pi*n - y - asin(y + 1.5) + %pi
endfunction

function plot_pereodical(f, x, n, pcolor, swap_axes)
    if ~exists("pcolor", "local") then pcolor = "b" end
    if ~exists("swap_axes", "local") then swap_axes = %F end
    for i=1:length(n)
        y = f(x, n(i))
        if swap_axes then
            plot(y, x, pcolor)
        else
           plot(x, y, pcolor)
        end
    end
endfunction


x1 = [-2.5:0.01:-0.5]'
x2 = [-4:0.01:4]'

plot(0.5 - cos(x2 - 0.5), x2, "k")
plot_pereodical(f1_p1, x1, -2:2, swap_axes=%T)
plot_pereodical(f1_p2, x1, -2:2, swap_axes=%T)
legend("x + cos(y - 0.5) - 0.5 = 0", "sin(x + y) - y - 1.5 = 0")

// Визуально можно определили начальное приближения для корня:
x01 = [1.4, -2.2]'

printf("------------ Поиск решения системы ------------")
[x11, x12, x13] = verbose_solver(x01, eps=10e-6)
norm(x11 - x12)
