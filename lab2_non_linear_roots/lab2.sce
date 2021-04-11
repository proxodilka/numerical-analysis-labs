format("e", 16);
function pretty_print(kwargs)
    to_print = ""
    for i=1:length(kwargs)
        sep = ""
        [key, value] = kwargs(i)(:)
        if value >= 0 then sep = " " end
        to_print = to_print + string(key) + "=" + sep + string(value) + " | "
    end
    printf(to_print + "\n")
endfunction

function [res]=f(x)
    // Исследуемая функция
    res = exp(-0.5*x) - 0.2 * x^2 + 1
endfunction

function [res]=df(x)
    // Производная исследуемой функции
    res = -0.5 * exp(-0.5*x) - 0.4 * x
endfunction

function [x0]=simple_iter(f, s, eps, x0, max_iter, verbose, check_convergence_cond, error_func, nprevs, verbose_func)
    // Находит корень функции f методом простых итераций.
    //
    // Параметры
    // ---------
    // f: функция,
    //    Функция, у которой нужно найти корень
    // s: функция (опционально, поумолчанию: s(x) = x + f(x)),
    //    Функция, принимающая на вход x^k и возвращающая x^(k+1) - более
    //    приближенное значение x*.
    // eps: скаляр (опционально, поумолчанию eps = 0.01),
    //    Требуемая точность корня.
    // x0: скаляр (опционально, поумолчанию x0 = 1),
    //    Начальное приближение корня.
    // max_iter: скаляр (опционально, поумолчанию max_iter = 1000),
    //    Максимальное число итераций метода.
    // verbose: булево (опционально, поумолчанию verbose = %F),
    //    Печатать шаги метода или нет
    // check_convergence_cond: булево, поумолчанию = %F,
    //    Проверять ли достаточные условия сходимости.
    //    ВНИМАНИЕ: при проверке условий сходимости будет вычислена первая производная функции.
    // error_func: функция, поумолчанию = `error_func(x_prevs) = abs(x_prevs(1) - x_prevs(2))`
    //    Функция, считающая оценку погрешности между текущим найденным решением и истиным корнем.
    //    Будет использоваться в условии остановки как: `if error_func(x_prevs) < eps then break`
    // nprevs: целое число, поумолчанию 1
    //    Кол-во предыдущих шагов, которые нужно передать в функцию `s(x)`. В функцию `error_func` 
    //    передается `nprevs + 1` шагов. Например, при `nprevs=2` в s(x) будут передаваться: `x^k, x^(k-1)`,
    //    а в `error_func`: `x^(k+1), x^k, x^(k-1)`.
    // verbose_func: функция, поумолчанию = `verbose_func(i, x, err) = printf(..., i, x, err)`
    //    Функция, принимающая номер текущего шага, текущее решение и текущую ошибку, и выводящая информацию
    //    об этом шаге на экран.
    //
    // Возвращает
    // ----------
    // Скаляр,
    //    Найденное значение x', удовлетворяющее заданной точности.
    //
    // Примечание
    // ----------
    // Метод простой итерации уточняет корень на каждом своем шаге, используя функцию s(x).
    // Более формально строится ряд:
    //     x^(k+1) = s(x^k) | s(x) = x + p(x) * f(x)
    // который при стремлении k->inf, должен сходится к точному решению.
    
    if ~exists("s", "local") then
        // Определим функцию s(x) поумолчанию: в качестве ф-ии p(x) возьмём константу 1
        function [_res]=s(x)
            xk = x(1)
            _res = xk + f(xk) 
        endfunction
    end
    if ~exists("x0", "local") then x0 = 1 end
    if ~exists("max_iter", "local") then max_iter = 1000 end
    if ~exists("eps", "local") then eps = 0.01 end
    if ~exists("verbose", "local") then verbose = %F end
    if ~exists("check_convergence_cond", "local") then check_convergence_cond = %F end
    if ~exists("nprevs", "local") then nprevs = 1 end
    if ~exists("error_func", "local") then 
        function [_res]=error_func(x)
            [x1, x0] = x(1:2)
            _res = abs(x0 - x1) 
        endfunction 
    end
    if ~exists("verbose_func", "local") then
        function verbose_func(i, x, err)
            printf("Шаг %i: ", i)
            pretty_print(list(list("x", x), list("f(x)", f(x)), list("error", err)))
        endfunction
    end

    if check_convergence_cond then
        // Проверим достаточные условия сходимости:
        //     |s(x)'| < 1  | x ∈ [a, b]
        // если это так, то x^0 = a и метод гарантированно сходится к x*
        if max(abs(numderivative(s, x0:0.1:x0 + 1))) >= 1 then
            printf("WARNING: Не выполнены достаточные условия сходимости метода!\n")
        end
    end

    x_prevs = list()
    // Заполняем массив с предыдущими значениями X значением поумолчанию
    for i=1:nprevs
        x_prevs(i) = %nan
    end
    
    for i=1:1:max_iter
        for j=2:nprevs
            x_prevs(j) = x_prevs(j - 1)
        end
        x_prevs(1) = x0
        x0 = s(x_prevs)
        err = error_func(list(x0, x_prevs(:)))
        if verbose 
            then verbose_func(i, x0, err)
        end
        if err < eps
            then break
    end
end
endfunction

function [res]=newtoon_method(f, df, eps, x0, max_iter, verbose, check_convergence_cond, root_locale)
    // Находит корень функции f методом Ньютона.
    //
    // Параметры
    // ---------
    // f: функция,
    //    Функция, у которой нужно найти корень
    // df: функция (опционально),
    //    Производная функции f. Если не переданна, производная будет
    //    вычислена автоматически.
    // eps: скаляр (опционально, поумолчанию eps = 0.01),
    //    Требуемая точность корня.
    // x0: скаляр (опционально, поумолчанию x0 = 1),
    //    Начальное приближение корня.
    // max_iter: скаляр (опционально, поумолчанию max_iter = 1000),
    //    Максимальное число итераций метода.
    // verbose: булево (опционально, поумолчанию verbose = %F),
    //    Печатать шаги метода или нет
    // check_convergence_cond: булево, поумолчанию = %F,
    //    Проверять ли достаточные условия сходимости, для этого дополнительно
    //    нужно передать параметр `root_locale`.
    //    ВНИМАНИЕ: при проверке условий сходимости будет вычислена вторая производная функции.
    // root_locale: list, опционально
    //    Список состоящий из двух элементов - начало и конец отрезка, на котором локализован
    //    корень. Используется только для проверки условий сходимости.
    //
    // Возвращает
    // ----------
    // Скаляр,
    //    Найденное значение x', удовлетворяющее заданной точности.
    //
    // Примечание
    // ----------
    // Метод Ньютона это итерационный метод, который на каждой итерации берет за
    // новое уточненное значение корня точку пересечения касательной к искомой функции
    // с осью OX.
    
    if ~exists("df", "local") then
        function [_res]=default_df(x)
            _res = diag(numderivative(f, x))
        endfunction
        df = default_df
    end
    if ~exists("x0", "local") then x0 = 1 end
    if ~exists("max_iter", "local") then max_iter = 1000 end
    if ~exists("eps", "local") then eps = 0.01 end
    if ~exists("verbose", "local") then verbose = %F end

    if check_convergence_cond then
        approx_domain = root_locale(1):0.1:root_locale(2)
        ddf_val = diag(numderivative(df, approx_domain))
        max_ddf = max(ddf_val)
        min_ddf = min(ddf_val)
        condition = (min_ddf * max_ddf >= 0 & (sum(isinf([max_ddf, min_ddf])) + sum(isnan([max_ddf, min_ddf]))) == 0 & prod(df(approx_domain)) ~= 0)
        if ~condition then
            printf("WARNING: Не выполнены достаточные условия сходимости метода!\n")
        end
    end

    function [_res]=s(x)
        xk = x(1)
        _res = xk - (f(xk) / df(xk))(:, 1)
    endfunction
    
    function [_res]=_error_func(x_prevs)
        [x2, x1, x0] = x_prevs(1:3)
        qn = (x2 - x1) / (x1 - x0)
        _res = abs((x2 - x1) / (1 - qn))
    endfunction
    
    function _verbose_func(i, x, err)
        printf("Шаг %i: ", i)
        pretty_print(list(list("x", x), list("f(x)", f(x)), list("f''(x)", df(x)), list("error", err)))
    endfunction
    
    res = simple_iter(f, s, eps, x0, max_iter, verbose, error_func=_error_func, verbose_func=_verbose_func, nprevs=2, check_convergence_cond=%F)
endfunction

function [res]=secant_method(f, eps, x0, x1, max_iter, verbose)
    if ~exists("eps", "local") then eps = 0.01 end
    if ~exists("x0", "local") then x0 = 1 end
    if ~exists("x1", "local") then x1 = x0 + eps end
    if ~exists("max_iter", "local") then max_iter = 1000 end
    if ~exists("verbose", "local") then verbose = %F end
    
    function [_res]=s(x)
        [_x1, _x0] = x(1:2)
        if isnan(_x0) then
            _x0 = x0
            _x1 = x1
        end
        _res = _x1 - ((_x1 - _x0) / (f(_x1) - f(_x0))) * f(_x1)
    endfunction
    
    res = simple_iter(f, s, eps, x0, max_iter, verbose, nprevs=2)
endfunction

// Изобразим функцию на графике
subplot(121)
acc = 0.01
x = [-5:acc:5]'
plot2d(x, f(x))

// Разобьем исходную функцию на две, их тоже нарисуем.
function [res]=f1(x)
    res = exp(-0.5 * x)
endfunction

function [res]=f2(x)
    res = 0.2 * x ^ 2 - 1
endfunction

subplot(122)
plot2d(x, [f1(x) f2(x)], [1, 2])

// Из графиков видно, что как минимум один корень располагается на отрезке [2, 3]
// Из производной исследуемой функции также видно, что её знак всегда отрицателен
//     f(x)' = -0.5*e^(-0.5x) - 0.4x < 0  (при x > 0: очевидно, при x < 0: значение экспоненты
//                                         по модулю, всегда будет больше значения линейной функции)
// => можно заключить что на всей числовой прямой у функции будет всего один корень, да и тот
// находится на отрезке [2, 3]

printf("\n\n------------ Метод простых итераций ------------\n")
x = simple_iter(f, eps=0.001, x0=2, verbose=%T, check_convergence_cond=%T)
printf("\nРешение, полученное методом простых итераций:\n\tx = %.8f | f(x) = %.8f", x, f(x))

printf("\n\n------------ Метод Ньютона ------------\n")
x = newtoon_method(f, x0=2, eps=1.e-6, verbose=%T, check_convergence_cond=%T, root_locale=[2, 4])
printf("\nРешение, полученное методом Ньютона:\n\tx = %.8f | f(x) = %.8f", x, f(x))

printf("\n\n------------ Метод секущих ------------\n")
x = secant_method(f, eps=1.e-6, x0=2, x1=3, verbose=%T)
printf("\nРешение, полученное методом секущих:\n\tx = %.8f | f(x) = %.8f", x, f(x))

printf("\n\n------------ fsolve ------------\n")
x = fsolve([1], f)
printf("\nРешение, полученное методом fsolve:\n\tx = %.8f | f(x) = %.8f", x, f(x))





