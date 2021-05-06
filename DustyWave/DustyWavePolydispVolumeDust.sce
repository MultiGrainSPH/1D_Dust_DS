//НАЧАЛЬНЫЕ ДАННЫЕ
rho_g = 1; //начальная плотность газа
N = 2; //количество фракций пыли
rho_i = [1, 1]; //массив начальной плотности пыли//все rho_i больше 0 и их сумма меньше rho_s
t_i = [0.01, 0.001]; //массив времён релаксации
rho_s = 4; //истинная плотность пыли
Cs = 1.0; //скорость звука в газе
k = 2*%pi; //волновое число

delta_rho_g_Wave = 1;
A_V=10^(-4)*Cs; //амплитуда для скорости
A_D=10^(-4)*rho_g; //амплитуда для плотности

С_Const = 1;//длина отрезка координат
number_Of_Points = 100;//количество узлов, по которым строим графики
t_Selected = 0.1;//в какой момент времени рассматриваем волну (t*)

path='D:\Program_for_Article\Output files';//путь до текстовых файлов

color_Arr=['blue', 'green'];
//КОНЕЦ НАЧАЛЬНЫХ ДАННЫХ

//
theta_i = rho_i / rho_s;
Theta = sum(theta_i);
//disp(theta_i);
//disp(Theta);
for i=1:N
    rho_div_t(i) = rho_i(i) / t_i(i);
end

//ЗАДАЁМ МАТРИЦУ (7)
for i=1:2*N+2
    for j=1:2*N+2
        M(i,j) = 0;
    end
end
//первая четверть (1,1)->(N+1,N+1)
M(1,1) = poly([sum(rho_div_t),-rho_g],'w','c');
for i=1:N
    M(1, i+1) = -rho_div_t(i);
    M(i+1, 1) = -rho_div_t(i);
    M(i+1,i+1) = poly([rho_div_t(i),-rho_i(i)],'w','c');
end
//вторая четверть (1,N+2)->(N+1,2N+1)
M(1,N+2) = Cs^2 * imult(k);
for i=1:N
    M(1,i+N+2) = Cs^2 * imult(k) * rho_g / (rho_s * (1 - Theta));
end
for i=1:N
    M(i+1,N+2) = (theta_i(i) / (1 - Theta)) * M(1,N+2);
    for j=1:N
        M(i+1,j+N+2) = (theta_i(i) / (1 - Theta)) * M(1,j+N+2);
    end
end
//третья четверть (N+2,1)->(2N+1,N+1)
M(N+2,1) = imult(k);
for i=1:N
    M(i+N+2,i+1) = imult(k);
end
//четвёртая четверть (N+2,N+2)->(2N+2,2N+2)
M(N+2,N+2) = poly([0,-rho_g^(-1)],'w','c');
for i=1:N
    M(i+N+2,i+N+2) = poly([0,-rho_i(i)^(-1)],'w','c');
end

//disp(det(M));

for i=1:2*N+2
    arr_Roots(i) = roots(det(M))(i);
end
disp(arr_Roots);

//ВЫБИРАЕМ КОРЕНЬ w_root ПО НАИБОЛШЕЙ МНИМОЙ ЧАСТИ
w_root = arr_Roots(1);
for i=2:2*N+2
    if imag(arr_Roots(i)) > imag(w_root) then
        w_root = arr_Roots(i);
    end
end
disp(w_root);

//ЗАДАЁМ МАТРИЦУ (7) БЕЗ СТОЛБЦА И СТРОКИ
for i=1:2*N+1
    for j=1:2*N+1
        new_M(i,j) = 0;
    end
end
//первая четверть
for i=1:N
    new_M(i, 1) = -rho_div_t(i);
    new_M(i,i+1) = rho_div_t(i) - (rho_i(i) * w_root);
end
//вторая четверть //
for i=1:N
    for j=1:N
        new_M(i,j+N+1) = (theta_i(i) / (1 - Theta)) * Cs^2 * imult(k) * rho_g / (rho_s * (1 - Theta));
    end
end
//третья четверть
new_M(N+1,1) = imult(k);
for i=1:N
    new_M(i+N+1,i+1) = imult(k);
end
//четвёртая четверть
for i=1:N
    new_M(i+N+1,i+N+1) = -w_root / rho_i(i);
end

//ЗАДАЁМ ПРАВУЮ ЧАСТЬ
for i=1:N
    B(i) = (theta_i(i) / (1 - Theta)) * Cs^2 * imult(k) * (-delta_rho_g_Wave);
end
B(N+1) = -w_root / rho_g * (-delta_rho_g_Wave);
for i=1:N
    B(i+N+1) = 0;
end

//НАХОДИМ РЕШЕНИЕ СЛАУ
X = inv(new_M) * B; //массив с дельтами скоростей и плотностей пыли
//disp(X);

//ЗАПИСЫВАЕМ ЗНАЧЕНИЯ ВОЗМУЩЕНИЙ С ВОЛНАМИ
delta_v_Wave = X(1);
for i=1:N
    delta_ui_Wave(i) = X(i+1);
end
for i=1:N
    delta_rhoi_Wave(i) = X(i+N+1);
end

//ВОЗМУЩЕНИЯ
function y=real_delta_f(delta_f_Wave, x, t, A)
    y = A * delta_f_Wave * exp(imult(k) * x - w_root * t);
    y = real(y);
endfunction

//СЕТКА
function X=uniform_Grid(a,b,N)//равномерная сетка на [a,b] с N узлами
    h=(b-a)/(N-1);
    for n=1:N
        X(n)=a+h*(n-1);
    end 
endfunction

X_x=uniform_Grid(0,С_Const,number_Of_Points);

//СОЗДАЁМ ИМЕНА ФАЙЛАМ ДЛЯ ГАЗА И ПЫЛИ ПРИ t=0, t=t*
for i=1:1+N
    if i==1 then
        name_Arr(i) = path + '\' + 'gas.txt';
    else
        name_Arr(i) = path + '\' + 'dust ' + string(i-1) + '.txt';
    end
end

//СОЗДАЁМ ПУСТЫЕ ФАЙЛЫ ДЛЯ ГАЗА И ПЫЛИ И ЗАПИСЫВАЕМ В НИХ ШАПКУ ТАБЛИЦЫ
for i=1:1+N
    f_w = mopen(name_Arr(i), 'wt');
    mfprintf(f_w, '#    x    Density (t = 0)    Density (t = ' + string(t_Selected) + ')    Velocity (t = 0)    Velocity (t = ' + string(t_Selected) + ')\n');
    mclose(f_w);
end

//ЗАПИСЫВАЕМ ДАННЫЕ В ФАЙЛЫ
f_a = mopen(name_Arr(1), 'at');
for i=1:length(X_x)
    Y1 = rho_g + real_delta_f(delta_rho_g_Wave, X_x(i), 0, A_D);
    Y2 = rho_g + real_delta_f(delta_rho_g_Wave, X_x(i), t_Selected, A_D);
    Y3 = real_delta_f(delta_v_Wave, X_x(i), 0, A_V);
    Y4 = real_delta_f(delta_v_Wave, X_x(i), t_Selected, A_V);
    mfprintf(f_a, '%f %f %f %f %f\n', X_x(i), Y1, Y2, Y3, Y4);
end
mclose(f_a);

for s=1:N
    f_a = mopen(name_Arr(s+1), 'at');
    for i=1:length(X_x)
        Y1 = rho_i(s) + real_delta_f(delta_rhoi_Wave(s), X_x(i), 0, A_D);
        Y2 = rho_i(s) + real_delta_f(delta_rhoi_Wave(s), X_x(i), t_Selected, A_D);
        Y3 = real_delta_f(delta_ui_Wave(s), X_x(i), 0, A_V);
        Y4 = real_delta_f(delta_ui_Wave(s), X_x(i), t_Selected, A_V);
        mfprintf(f_a, '%f %f %f %f %f\n', X_x(i), Y1, Y2, Y3, Y4);
    end
    mclose(f_a);
end


//НАХОДИМ ЗНАЧЕНИЯ ПЛОТНОСТИ И СКОРОСТИ ГАЗА В УЗЛАХ
for i=1:length(X_x)
    Y_x_rho_g_0(i) = rho_g + real_delta_f(delta_rho_g_Wave, X_x(i), 0, A_D);
    Y_x_rho_g(i) = rho_g + real_delta_f(delta_rho_g_Wave, X_x(i), t_Selected, A_D);
    Y_x_v_g_0(i) = real_delta_f(delta_v_Wave, X_x(i), 0, A_V);
    Y_x_v_g(i) = real_delta_f(delta_v_Wave, X_x(i), t_Selected, A_V);
end

//РИСУЕМ ГРАФИКИ
subplot(1,2,1);
plot2d(X_x, Y_x_rho_g_0, style=color(color_Arr(1)), strf='181');
plot2d(X_x, Y_x_rho_g, style=color(color_Arr(2)), strf='181');
hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
g=get('current_axes');
g.title.text='Gas Density';

subplot(1,2,2);
plot2d(X_x, Y_x_v_g_0, style=color(color_Arr(1)), strf='181');
plot2d(X_x, Y_x_v_g, style=color(color_Arr(2)), strf='181');
hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
g=get('current_axes');
g.title.text='Gas Velocity';
