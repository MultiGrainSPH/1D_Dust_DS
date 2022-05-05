//НАЧАЛЬНЫЕ ДАННЫЕ
rho_g = 1; //начальная плотность газа
e_g=1;
N = 1; //количество фракций пыли
rho_i = [1]; //массив начальной плотности пыли//все rho_i больше 0 и их сумма меньше rho_s
e_i=[1]
t_i = [10]; //массив времён релаксации
rho_s = 4; //истинная плотность пыли
Cs = 2/3; //скорость звука в газе
k = 2*%pi; //волновое число
g=4/3;
Cp=1.0;
Cv=1.0;
Cdv=1.0;
zt=10;


delta_rho_g_Wave = 1;
A_V=Cs; //амплитуда для скорости
A_D=rho_g; //амплитуда для плотности

С_Const = 1;//длина отрезка координат
number_Of_Points = 100;//количество узлов, по которым строим графики
t_Selected = 0.1;//в какой момент времени рассматриваем волну (t*)

path='C:\Program_for_Article\Output files';//путь до текстовых файлов

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
disp(rho_i(1))
//ЗАДАЁМ МАТРИЦУ (7)


for i=1:3*N+3
    for j=1:3*N+3
        M(i,j) = 0;
    end
end

M(1,1)=poly([0,imult(1)],'w','c');
M(1,2)=0;
M(1,3)=-k*rho_g;
M(1,4)=0;
M(1,5)=0;
M(1,6)=0;
M(2,1)=0;
M(2,2)=poly([0,imult(1)],'w','c');
M(2,3)=0;
M(2,4)=-k*rho_i(1);
M(2,5)=0;
M(2,6)=0;
M(3,1)=-k*(g-1)*e_g;
M(3,2)=0;
M(3,3)=poly([imult(rho_div_t(1)),rho_g*imult(1)],'w','c');
M(3,4)=-imult(rho_div_t(1));
M(3,5)=-k*(g-1)*rho_g;
M(3,6)=0;
M(4,1)=0;
M(4,2)=0;
M(4,3)=-imult(rho_div_t(1));
M(4,4)=poly([imult(rho_div_t(1)),rho_i(1)*imult(1)],'w','c');
M(4,5)=0;
M(4,6)=0;
M(5,1)=0;
M(5,2)=0;
M(5,3)=-k*(g-1)*e_g*rho_g;
M(5,4)=0;
M(5,5)=poly([imult(rho_i(1)*Cp/Cv/zt),rho_g*imult(1)],'w','c');
M(5,6)=-imult(rho_i(1)*Cp/Cdv/zt);
M(6,1)=0;
M(6,2)=0;
M(6,3)=0;
M(6,4)=0;
M(6,5)=-imult(rho_i(1)*Cp/Cv/zt);
M(6,6)=poly([imult(rho_i(1)*Cp/Cdv/zt),rho_i(1)*imult(1)],'w','c');
disp(M);
disp(det(M));
for i=1:3*N+3
    arr_Roots(i) = roots(det(M))(i);
end
disp(arr_Roots);

//ВЫБИРАЕМ КОРЕНЬ w_root ПО НАИБОЛШЕЙ МНИМОЙ ЧАСТИ
w_root = arr_Roots(1);
for i=2:3*N+3
    if imag(arr_Roots(i)) > imag(w_root) then
        w_root = arr_Roots(i);
    end
end
disp(w_root);

//ЗАДАЁМ МАТРИЦУ (7) БЕЗ СТОЛБЦА И СТРОКИ
for i=1:3*N+2
    for j=1:3*N+2
        new_M(i,j) = 0;
    end
end


 
new_M(1,1)=0;
new_M(1,2)=-k*rho_g;;
new_M(1,3)=0;
new_M(1,4)=0;
new_M(1,5)=0;
new_M(2,1)=imult(w_root);
new_M(2,2)=0;
new_M(2,3)=-k*rho_i(1);
new_M(2,4)=0;
new_M(2,5)=0;
new_M(3,1)=0;
new_M(3,2)=imult(rho_div_t(1))+rho_g*imult(1)*w_root;
new_M(3,3)=-imult(rho_div_t(1));
new_M(3,4)=-k*(g-1)*rho_g;
new_M(3,5)=0;
new_M(4,1)=0;
new_M(4,2)=-imult(rho_div_t(1));
new_M(4,3)=imult(rho_div_t(1))+rho_i(1)*imult(1)*w_root;
new_M(4,4)=0;
new_M(4,5)=0;
new_M(5,1)=0;
new_M(5,2)=0;
new_M(5,3)=0;
new_M(5,4)=imult(rho_i(1)*Cp/Cv/zt+rho_g*w_root);
new_M(5,5)=-imult(rho_i(1)*Cp/Cdv/zt);
disp(new_M);
disp(inv(new_M));

//ЗАДАЁМ ПРАВУЮ ЧАСТЬ
B(1)=imult(1)*w_root;
B(2)=0;
B(3)=-k*(g-1)*e_g;
B(4)=0;
B(5)=0;
disp(B);

//НАХОДИМ РЕШЕНИЕ СЛАУ
 X = inv(new_M) * B;//массив с дельтами скоростей и плотностей пыли
disp(X);
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
        Y1 =  rho_i(s)+real_delta_f(delta_rhoi_Wave(s), X_x(i), 0,A_D);
        Y2 =  rho_i(s)+real_delta_f(delta_rhoi_Wave(s), X_x(i), t_Selected, A_D);
        Y3 = real_delta_f(delta_ui_Wave(s), X_x(i), 0, A_V);
        Y4 = real_delta_f(delta_ui_Wave(s), X_x(i), t_Selected, A_V);
        mfprintf(f_a, '%f %f %f %f %f\n', X_x(i), Y1, Y2, Y3, Y4);
    end
    mclose(f_a);
end


//НАХОДИМ ЗНАЧЕНИЯ ПЛОТНОСТИ И СКОРОСТИ ГАЗА В УЗЛАХ
for i=1:length(X_x)
    Y_x_rho_g_0(i) = real_delta_f(delta_rho_g_Wave, X_x(i), 0, A_D);
    Y_x_rho_g(i) = real_delta_f(delta_rho_g_Wave, X_x(i), t_Selected, A_D);
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

