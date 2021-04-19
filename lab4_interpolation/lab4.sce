x = [0.079, 0.637, 1.345, 2.095, 2.782]'
y = [-4.308, -0.739, 1.697, 4.208, 6.203]'

x_unknown = [x(1) + x(2), (x(1) + x(2)) / 2, (x(4) + x(5)) / 2, (x(1) + x(5)) / 3]'
x_approx = linspace(0, 3, n=100)'

function [res]=lagrange_poly(x, src_x, src_y)
    function [_res]=build_node(i, _x)
       [top_accum, bottom_accum] = list(1, 1)(:)
        for j=1:length(src_x)
            if j ~= i then
                top_accum = top_accum .* (_x - ones(_x) * src_x(j)) 
                bottom_accum = bottom_accum * (src_x(i) - src_x(j))
            end
        end
        _res = top_accum / bottom_accum
    endfunction
      
    res = zeros(x)
    for i=1:length(src_x)
        res = res + build_node(i, x) * src_y(i)
    end
endfunction

function [res]=nan_table(ndim)
    res = []
    for i=1:ndim
        for j=1:ndim
            res(i, j) = %nan
        end
    end
endfunction

function [res]=merge_tables(tables)
    res = nan_table(size(tables(1))(1))
    for i=1:length(tables(1))
        for j=1:length(tables)
            if ~isnan(tables(j)(i)) then res(i) = tables(j)(i) end
        end
    end
endfunction

function [res]=recursive_diff(x, y, as_table, shape)
    if ~exists("as_table") then as_table = %F end
    if ~exists("shape") then shape = 1 end
    if type(as_table) == type(%T) && as_table == %T then 
        as_table = list(length(x), length(x))
        shape = length(x)
    end
    if length(x) == 1 then
        res = y(1)
        if type(as_table) == type(list()) then
            res = nan_table(shape)
            res(as_table(:)) = y(1) 
        end
        return
    end
    
    if type(as_table) == type(list()) then
        l_idx = list(as_table(1), as_table(2) - 1)
        r_idx = list(as_table(1) - 1, as_table(2) - 1)

        l_diff = recursive_diff(x(2:$), y(2:$), l_idx, shape)
        r_diff = recursive_diff(x(1:$-1), y(1:$-1), r_idx, shape)
        res = merge_tables(list(r_diff, l_diff))

        l_value = res(l_idx(:))
        r_value = res(r_idx(:))
        res_value = (r_value - l_value) / (x($) - x(1))
        
        res(as_table(:)) = res_value
    else
        l_diff = recursive_diff(x(2:$), y(2:$))
        r_diff = recursive_diff(x(1:$-1), y(1:$-1))
        
        res = (l_diff - r_diff) / (x($) - x(1))
    end
endfunction

function [res]=newtoon_poly(x, x_src, y_src)
    res = 0
    for i=1:length(x_src)
        node_accum = 1
        for j=1:i - 1
            node_accum = node_accum .* (x - ones(x) * x_src(j))
        end
        res = res + node_accum .* recursive_diff(x_src(1:i), y_src(1:i))
    end
endfunction

function [res]=append_to_head(val, arr)
    for i=length(arr):-1:1
        arr(i + 1) = arr(i)
    end
    arr(1) = val
    res = arr
endfunction

function [res]=err(x, X, Y)
    accum = 1
    for i=1:length(X) - 1
        accum = accum * (x - X(i))
    end
    res = recursive_diff(X, Y) * accum
endfunction

lagr_interpl = lagrange_poly(x_approx, x, y)
lagr_predicted = lagrange_poly(x_unknown, x, y)
disp("Значения полинома Лагранжа в точках:", [
    "x(1) + x(2): ", "(x(1) + x(2)) / 2: ", "(x(4) + x(5)) / 2: ", "(x(1) + x(5)) / 3: ";
    string(lagr_predicted')
])

disp("Таблица раздельных расностей:", recursive_diff(x, y, as_table=%T))
newtoon_interpl = newtoon_poly(x_approx, x, y)

disp("Оценка погрешности интерполяционного полинома в точке (x(1) + x(2)) / 2:", err((x(1) + x(2)) / 2, x, y))

linear_interpl = interpln([x'; y'], x_approx)
spline_interpl = interp(x_approx, x, y, splin(x, y))

// Первый график: Лагран + Ньютон
scatter(x, y)
plot(x_approx, [lagr_interpl], 'LineWidth', 4)
plot(x_approx, [newtoon_interpl], '--r', 'LineWidth', 4)
legend("известные значения функции", "полином Лагранжа", "полином Ньютона", pos=4)
scf(1)

// Второй график: линейный + кубический сплайн
scatter(x, y)
plot(x_approx, [linear_interpl], "k", "LineWidth", 2)
plot(x_approx, [spline_interpl], "c", "LineWidth", 2)
legend("известные значения функции", "линейный сплайн", "кубический сплайн", pos=4)
scf(2)

// Третий график: все методы (Лагранж, Ньютон, Линейный сплайн, Кубический сплайн)
scatter(x, y)
plot(x_approx, [lagr_interpl], 'LineWidth', 4)
plot(x_approx, [newtoon_interpl], '--r', 'LineWidth', 4)
plot(x_approx, [linear_interpl], "k", "LineWidth", 2)
plot(x_approx, [spline_interpl], "c", "LineWidth", 2)
legend("известные значения функции", "полином Лагранжа", "полином Ньютона", "линейный сплайн", "кубический сплайн", pos=4)


