X=[0.351, 0.664, 0.978, 1.291, 1.605, 1.918, 2.232, 2.546, 2.859]
Y=[0.605, 0.265, 0.064, 0.116, 0.415, 0.728, 1.673, 3.138, 5.092];
x_approx = 0:0.1:3

function [y]=square_sum(x)
    y = 0
    for i=1:length(x)
        y = y + x(i) ^ 2
    end
endfunction

function [y]=approx(params, x)
    a = params(1)
    b = params(2)
    c = params(3)
    y = approx_fn(x, a, b, c)
endfunction

function [err]=error_fn(c, params)
    x = params(1, :)
    y = params(2, :)
    err = (y - approx(c, x))'
endfunction

function [y]=f1(x, a, b, c)
    y = a * x^2 + b*x + c
endfunction

function [y]=f2(x, a, b, c)
    y = a * x^(-2) + b * x^(-1) + c
endfunction

function [y]=f3(x, a, b, c)
    y = b * x ^ a + c
endfunction

function [y]=f4(x, a, b, c)
    y = b * exp(a * x) + c
endfunction

function [y]=f5(x, a, b, c)
    y = b * ((x + a)^-1) + c
endfunction

function [y]=f6(x, a, b, c)
    y = a * x + b * exp(-x) + c
endfunction

function [y]=f7(x, a, b, c)
    y = a * (x ^ -1) + b * exp(x) + c
endfunction

function [y]=f8(x, a, b, c)
    y = a * log(x) + b * exp(x) + c
endfunction

function [y]=f9(x, a, b, c)
    y = b * exp(-a * (x + c) ^ 2) + c
endfunction

function [y]=f10(x, a, b, c)
    y = a * x ^ (0.5) + b * sin(x) + c
endfunction

f=get("current_figure")
f.figure_size=[720,541]
fn_list = list(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10)
// fn_list = list(f9, f10)
// fn_list = list()
for i=1:length(fn_list)
    approx_fn = fn_list(i)
    [p, err] = datafit(error_fn, [X; Y], [1, 1, 0])
    printf("\n---------- Функция f%i ----------\n", i)
    disp("Найденные коэффициенты:", p)
    printf("\nОценка погрешности datafit: %e\n", err)
    predicted_vals = error_fn(p, [X; Y])
    disp("Значение невязок в исходных точках:", predicted_vals' * predicted_vals)
    disp("Средне-квадратичное отклонение ошибки:", stdev(predicted_vals))
    // scf(i)
    subplot(2, 5, i)
    scatter(X, Y)
    plot(x_approx, approx(p, x_approx))
    legend("известные значения", "аппроксимация функцией f" + string(i))
end

function [A, b]=build_system(bias_fns, x, y)
    A = eye(length(bias_fns), length(bias_fns))
    b = eye(length(bias_fns))
    
    for i=1:length(bias_fns)
        for j=1:length(bias_fn)
            // disp(bias_fns(j)(x))
            A(i, j) = bias_fns(i)(x) * (bias_fns(j)(x))'
        end
        b(i) = bias_fns(i)(x) * y'
    end
endfunction
// Определим базисные функции для f8
printf("\n---------- Оптимальные коэффициенты функция f8 ----------\n", i)
bias_fn = list(log, exp, ones)
approx_fn = fn_list(8)
[A, b] = build_system(bias_fn, X, Y)
p = linsolve(A, -b)
disp("Оптимальные параметры для модели 8:", p)
predicted_vals = error_fn(p, [X; Y])
disp("Значение невязок в исходных точках:", predicted_vals' * predicted_vals)
disp("Средне-квадратичное отклонение ошибки:", stdev(predicted_vals))





