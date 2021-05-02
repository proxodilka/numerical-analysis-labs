// format("e", 14)
function [x, y]=transform(A, b)
    order = []
    for row=1:1:size(A)(1)
        [el, k] = max(A(row, :))
        order(row) = k
    end 
    x = A(order, :)
    y = b(order, :)
endfunction

function [rel_eps, abs_eps]=get_eps(ref_x, res_x)
    abs_eps = norm(res_x - ref_x)
    rel_eps = abs_eps / norm(ref_x)
endfunction

function [rel_err, abs_err]=get_err(A, b, err)
    cnd = cond(A)
    diff_A = err / norm(A)
    diff_B = err / norm(b)
    
    rel_err = (cnd / (1 - cnd * diff_A)) * (diff_B + cnd * diff_A)
    abs_err = norm(A^-1) * err
endfunction

function [x]=verbose_solver(A, b, solver, method_name, args, ref_x)
    if ~exists("args", "local") then args = list() end
    if ~exists("ref_x") then ref_x = A^-1 * b end
    
    printf("\n\n------------ " + method_name + " ------------\n")
    x = solver(A, b, args(:), verbose=%T)
    disp(method_name + " выдал следующее решение:", x)
    err = norm(A * x - b)
    printf("Норма вектора невязки: %f\n", err)
    [rel_err, abs_err] = get_eps(ref_x, x)
    printf("Относительная погрешность решения: %.4e\n", rel_err)
    printf("Абсолютная погрешность решения: %.4e\n", abs_err)
endfunction

function [x]=solve_with_gaus(A, b, verbose)
    /// Решение методом Гауса:
    /// 1. Приводим систему к треугольному виду с помощью оператора `rref`
    /// 2. Последовательно выражаем неизвестные переменные
    triang = rref([A, b])

    triang_A = triang(:, 1:$-1)
    triang_B = triang(:, $)
    
    gaus_x = zeros(b)
    
    for i=size(b)(1):-1:1
        rsum = 0
        for j=i:size(b)(1)
           rsum = rsum + triang_A(i, j) * gaus_x(j)
        end
        gaus_x(i) = (triang_B(i) - rsum) / triang_A(i, i)
    end
    x = gaus_x
endfunction

function [x]=solve_with_iter(A, b, eps, max_iter, verbose)
    /// Решение методом простых итераций:
    /// 1. Приводим систему к виду, чтобы на главной диагонали были макс. элементы матрицы `A`
    /// 2. Приводим систему к виду: x = Q*x + c, где:
    ///     a) Q = E - D*A
    ///     b) c = D*b
    ///     c) max(spec(Q)) < 1
    ///    Тогда, рассмотрев рекурентный ряд:
    ///         x^(k+1) = Q*x^(k) + c
    ///    обнаружим, что он сходится и имеет пределом матрицу:
    ///         x' = (E - (E - D*A))^(-1) * c
    ///         A*x' = b  => x' - решение исходной системы
    if ~exists("eps", "local") then eps = 0.01 end
    if ~exists("max_iter", "local") then max_iter = 1000 end
    if ~exists("verbose", "local") then verbose = %F end
    
    [A, b] = transform(A, b)
    
    /// Построим новые матрицы `D*A` и `с` следующим образом:
    ///     D*A(i, j) = A(i, j)/A(i, i)
    ///     c(i) = b(i)/A(i, i)
    diags = diag(A)
    for i=1:1:size(A)(1)
        DA(i, :) = A(i, :) / diags(i)
        c(i) = b(i) / diags(i)
    end
    
    Q = (eye(DA) - DA)
    
    /// Проверяем НИД условие сходимости метода
    max_spec = max(abs(spec(Q)))
    if (max_spec >= 1) then
        printf("WARNING: Максимальное по модулю собственное число преобразованной матрицы не удовлетворяет НИД условию сходимости: %f >= 1", max_spec)
    end
    
    x = c
    /// Априорной оценкой точности решения является:
    ///      |x' - x^k| <= |Q^k|*|(E - Q)^-1 * c - x^0|
    /// Соответственно, отсюда можно вывести необходимое число шагов `k` для получения
    /// решения с заданной точностью
    k = int(log(eps / (norm(x) + norm(c) / (1 - norm(Q)))) / log(norm(Q)) + 1)
    if k > max_iter then
        printf("WARNING: Оценочное число итераций больше максимального: %i > %i, будет сделано $max_iter = %i итераций.", k, max_iter, max_iter)
        k = max_iter
    end
    printf("Запуск итерационного метода на %i шагов...\n", k)
    for i=1:k
        x_prev = x
        x = Q*x + c
        if verbose then
            err = norm(x - x_prev)
            printf("\nШаг: %i | Погрешность: %.8f |\nТекущее решение:", i, err)
            disp(x)
        end
    end
endfunction
    
function [x]=solve_with_zeidel(A, b, eps, max_iter, verbose)
    /// Решение методом Зейделя:
    /// Решение метода сходно с методом простых итераций, отличие в этот раз
    /// будет заключатся лишь в том, что при обновлении решения будут учитываться
    /// уже посчитанные компоненты вектора на этой итерации
    if ~exists("eps", "local") then eps = 0.01 end
    if ~exists("max_iter", "local") then max_iter = 1000 end
    if ~exists("verbose", "local") then verbose = %F end
    
    [A, b] = transform(A, b)

    diags = diag(A)
    for i=1:1:size(A)(1)
        DA(i, :) = A(i, :) / diags(i)
        c(i) = b(i) / diags(i)
    end
    
    Q = (eye(DA) - DA)
    
    /// Проверяем НИД условие сходимости метода
    max_spec = norm(Q)
    if (max_spec >= 1) then
        printf("WARNING: Значение нормы преобразованной матрицы не удовлетворяет НИД условию сходимости: %f >= 1", max_spec)
    end
 
    x = c   
    eps = eps * (1 - norm(Q)) / norm(Q)
    for i=1:max_iter
        /// Сохраняем предыдущий результат для подсчета условия остановки
        x_prev = x
        for j=1:size(x)(1)
            x(j) = sum(Q(j, :) * x) + c(j)
        end
        err = norm(x_prev - x)
        if verbose then
            printf("\nШаг: %i | Погрешность: %f.8 |\nТекущее решение:", i, err)
            disp(x)
        end
        if err < eps
            then break
        end
    end
endfunction

A = [
    0.197 0.219 0.274 3.127
    0.186 0.275 2.987 0.316
    0.329 2.796 0.179 0.278
    2.389 0.273 0.126 0.418

]

b = [
    0.869
    0.529
    0.297
    0.144
]

[A, b] = transform(A, b)
disp("Преобразованная система с преобладанием диагональных элементов:", [A, b])

cond_n = cond(A)
printf("Число обусловнености равно: %f\n", cond_n)

if (cond_n > 100) 
    then printf("WARNING: Число обусловленности слишком большое, точность последующего решения под вопросом\n")
end

err = 0.001

[rel_err, abs_err] = get_err(A, b, err)
printf("Оценка относительной погрешности решения: %f\n", rel_err)
printf("Оценка абсолютной погрешности решения: %f\n", abs_err)

verbose_solver(A, b, solve_with_gaus, "Метод Гаусса")
verbose_solver(A, b, solve_with_iter, "Метод простых итераций", list(0.01))
verbose_solver(A, b, solve_with_iter, "Метод Якоби", list(0.001))
verbose_solver(A, b, solve_with_zeidel, "Метод Зейделя", list(0.001))
