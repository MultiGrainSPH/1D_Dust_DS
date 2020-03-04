//НАЧАЛЬНЫЕ ДАННЫЕ//INITIAL DATA
rho_g=1;//плотность газа//gas density
rho_Arr=[0.1 0.233333 0.366667 0.5];//массив плотностей фракций//fraction density array
ts_Arr=[0.1 0.215443 0.464159 1];//массив времён релаксации фракций//array of fraction relaxation times

Cs_Const=1;//скорость звука в газе//speed of sound in gas

С_Const=1;//длина отрезка координат//coordinate length
n_Const=1;//количество волн на длину отрезка//number of waves per segment length
T_Const=2;//длина отрезка времени//time length

t_Selected=1;//в какой момент времени рассматриваем волну (t*)//at what point in time we consider the wave (t*)
x_Selected=0;//в какой координате рассматриваем волну (x0)//in which coordinate we consider the wave (x0)

type_Of_Plot='t';//вдоль чего меняется график: x-координаты, t-время//along which the graph is changing: x-coordinates, t-time

overwrite_Files='on';//перезапись файлов при отрисовке графиков (on-да)//rewrite files when drawing graphs (on-yes)

number_Of_Points=100;//количество узлов, по которым строим графики//the number of nodes on which to build graphs
A_Velocity=10^(-4)*Cs_Const;//амплитуда для скорости//amplitude for velocity
A_Density=10^(-4)*rho_g;//амплитуда для плотности//amplitude for density

ed=10^(-10);//точность нахождения корней полинома//accuracy of finding the roots of the polynomial

name='Initial data';//имя файла, в который записываем данные для статьи//the name of the file in which we write the data for the article
path='C:\Users\User\Desktop\ГАЗОПЫЛЬ\ПРОГРАММА ДЛЯ ИНСТИТУТА\Output files';//путь до текстовых файлов//path to text files

//ЦВЕТА ГРАФИКОВ//GRAPHIC COLORS
color_Arr_Article=['blue', 'orange', 'red', 'green', 'purple'];//для графиков из статьи//for graphs from an article
color_Arr_Double=['green', 'blue'];//для двойных графиков//for double graphы

//КОНЕЦ НАЧАЛЬНЫХ ДАННЫХ//END OF INITIAL DATA

//СОЗДАЕМ ПУСТОЙ ТЕКСТОВЫЙ ФАЙЛ//CREATE AN EMPTY TEXT FILE
if overwrite_Files=='on' then
    full_Name=path+'\'+name+'.txt';
    f_w=mopen(full_Name, 'wt');
    mclose(f_w);
end

//СОЗДАЁМ ИМЕНА ФАЙЛАМ ДЛЯ ГАЗА И ПЫЛИ ПРИ t=0, t=t* И x=x0//CREATE NAMES FOR GAS AND DUST FILES AT t=0, t=t* AND x=x0
if overwrite_Files=='on' then
    for i=1:(1+length(rho_Arr))
        if i==1 then
            name_Arr_Double(i)=path+'\'+'gas (t=0).txt';
            name_Arr_Double(i+1+length(rho_Arr))=path+'\'+'gas (t='+string(t_Selected)+').txt';
            name_Arr_Normalized(i)=path+'\'+'gas (x='+string(x_Selected)+').txt';
        else
            name_Arr_Double(i)=path+'\'+'dust '+string(i-1)+' (t=0).txt';
            name_Arr_Double(i+1+length(rho_Arr))=path+'\'+'dust '+string(i-1)+' (t='+string(t_Selected)+').txt';
            name_Arr_Normalized(i)=path+'\'+'dust '+string(i-1)+' (x='+string(x_Selected)+').txt';
        end
    end
end

//СЧИТАЕМ КОНСТАНТУ k И МАССИВ epsilon//CALCULATE K CONSTANT AND ARRAY epsilon
//смотри после (49)//look after (49)
k_Const=2*%pi*n_Const/С_Const;
epsilon_Arr=rho_Arr/rho_g;

//ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ ПОЛИНОМА//AUXILIARY FUNCTIONS FOR POLINOMA
//смотри (45)//see (45)
function y=polynom(x)//значение полинома в точке x//polynomial value at point x
    ws=k_Const*Cs_Const;
    L=length(epsilon_Arr);
    if L==1 then
        S=1-x*ts_Arr(1);
        y=x^2*(S+epsilon_Arr(1))+ws^2*S;
    else
        p=1;
        s=0;
        for i=1:L
            p=p*(1-x*ts_Arr(i));
            pi=1;
            for j=1:L
                if j~=i then
                    pi=pi*(1-x*ts_Arr(j));
                end
            end
            s=s+epsilon_Arr(i)*pi;
        end
        y=x^2*(p+s)+ws^2*p;
    end
endfunction

function Y=dichotomy(a,b,f,e)//поиск нуля функции f с точностью e между a и b// finding a zero function f with e accuracy between a and b
    if f(a)==0 then
        m=a;
    else
        while abs(f(a)-f(b))>e & abs(a-b)>e
            m=(b+a)/2;
            if f(a)*f(m)<0 then
               b=m;
            elseif f(a)*f(m)>0 then
                a=m;
            else
                a=m;
                b=m;
            end
        end
        m=(b+a)/2;
    end
    Y(1)=m;
    Y(2)=f(m);
endfunction

function Y=decrease(X)//сортировка по убыванию//sorting by descending
    L=length(X);
    Y=X;
    if L>1 then
        for i=1:L-1
            for j=i+1:L
                if Y(j)>Y(i) then
                    temp=Y(i);
                    Y(i)=Y(j);
                    Y(j)=temp;
                end
            end
        end
    end
endfunction

function y=opposite_Sign(f,x)//ищем точку y такую, что f(x)*f(y)<0//we are looking for a point y such that f(x)*f(y)<0
    y=x;
    while f(x)*f(y)>=0
        y=y+abs(x);
    end
endfunction

//НАХОЖДЕНИЕ НУЛЕЙ ПОЛИНОМА//FINDING ZEROS OF POLYNOMIAL
ts_Arr_Sort=decrease(ts_Arr);
for i=1:length(ts_Arr_Sort)
    if i<length(ts_Arr_Sort) then
        D=dichotomy(1/ts_Arr_Sort(i),1/ts_Arr_Sort(i+1),polynom,ed);
    else
        if polynom(1/ts_Arr_Sort(i))==0 then
            D=[1/ts_Arr_Sort(i), polynom(1/ts_Arr_Sort(i))];
        else
            inf_Point=opposite_Sign(polynom,1/ts_Arr_Sort(i));
            D=dichotomy(1/ts_Arr_Sort(i),inf_Point,polynom,ed);
        end
    end
    roots_Arr(i)=D(1);
end

//КОЭФФИЦИЕНТЫ a0 И aN//COEFFICIENTS a0 AND aN
//смотри (D13), (D14)//see (D13), (D14)
a0=(-1)^length(ts_Arr)*(k_Const*Cs_Const)^2;
aN=0;
for i=1:length(ts_Arr)
    a0=a0/ts_Arr(i);
    aN=aN-(1+epsilon_Arr(i))/ts_Arr(i);
end

//КОЭФФИЦИЕНТЫ p И q//COEFFICIENTS p AND q
//смотри (D11), (D12)//see (D11), (D12)
p=aN;
q=(-1)^length(roots_Arr)*a0;
for i=1:length(roots_Arr)
    p=p+roots_Arr(i);
    q=q/roots_Arr(i);
end

//ПРЯМОЙ И ОБРАТНЫЙ ПЕРЕВОД ЗНАКА МНИМОЙ ЧАСТИ В СТРОКУ//DIRECT AND REVERSE TRANSFER OF THE SIGN OF THE IMAGINARY PART IN A LINE
function y=imag_Sign(x)
    if imag(x)>=0 then
        y='+';
    else
        y='-';
    end
endfunction

function y=imag_Sign_Rev(x)
    if imag(x)>=0 then
        y='-';
    else
        y='+';
    end
endfunction

//МНИМЫЙ КОРЕНЬ//IMAGINARY ROOT
w_Root=roots([1,p,q])(2);

//СЧИТАЕМ ДЕЛЬТЫ ПЛОТНОСТЕЙ И СКОРОСТЕЙ С ВОЛНАМИ//WE COUNT DENSITS OF DENSITY AND SPEED WITH WAVES
//смотри (46), (47), (48)//see (46), (47), (48)
delta_rho_g_Wave=1;
for i=1:length(rho_Arr)
    delta_rho_d_Wave(i)=(1/(1-w_Root*ts_Arr(i)))*(delta_rho_g_Wave/rho_g)*rho_Arr(i);//(48)
end
ws_Const=k_Const*Cs_Const;
delta_v_g_Wave=-%i*(w_Root/ws_Const)*(delta_rho_g_Wave/rho_g)*Cs_Const;//(46)
for i=1:length(rho_Arr)
    delta_v_d_Wave(i)=-%i*(w_Root/ws_Const)*(1/(1-w_Root*ts_Arr(i)))*(delta_rho_g_Wave/rho_g)*Cs_Const;//(47)
end

//СОЗДАЁМ ПУСТЫЕ ФАЙЛЫ ДЛЯ ГАЗА И ПЫЛИ И ЗАПИСЫВАЕМ В НИХ ШАПКУ ТАБЛИЦЫ//WE CREATE EMPTY FILES FOR GAS AND DUST AND WRITE IN THEM THE BEGINNING OF THE TABLE
if overwrite_Files=='on' then
    if type_Of_Plot=='x' then
        for i=1:(1+length(rho_Arr))*2
            f_w=mopen(name_Arr_Double(i), 'wt');
            if i==1 | i==1+1+length(rho_Arr) then
                D_Wave=delta_rho_g_Wave;
                V_Wave=delta_v_g_Wave;
            else
                if i<=1+length(rho_Arr) then
                    D_Wave=delta_rho_d_Wave(i-1);
                    V_Wave=delta_v_d_Wave(i-1);
                else
                    D_Wave=delta_rho_d_Wave(i-1-1-length(rho_Arr));
                    V_Wave=delta_v_d_Wave(i-1-1-length(rho_Arr));
                end
            end
            
            mfprintf(f_w, '#      x            Density           Velocity\n');
            mclose(f_w);
        end
    elseif type_Of_Plot=='t' then
        for i=1:(1+length(rho_Arr));
           f_w=mopen(name_Arr_Normalized(i), 'wt');
           mfprintf(f_w, '# Normalized fraction density = density delta / initial density\n');
           mfprintf(f_w, '# Normalized fraction velocity = velocity delta / Cs\n');
           mfprintf(f_w, '#\n');
           
           mfprintf(f_w, '#      t     Normalized Density     Normalized Velocity\n');
           mclose(f_w);
        end
    end
end

//ЗАПИСЫВАЕМ ТАБЛИЧНЫЕ ДАННЫЕ В ФАЙЛ//WRITE TABLE DATA IN FILE
if overwrite_Files=='on' then
    f_a=mopen(full_Name, 'at');
    mfprintf(f_a, 'Cs\n');
    mfprintf(f_a, '  %f\n', Cs_Const);
    mfprintf(f_a, '\n');
    
    mfprintf(f_a, 'Initial densities\n');
    mfprintf(f_a, '  g: %f\n', rho_g);
    for i=1:length(rho_Arr)
        mfprintf(f_a, '  %d: %f\n', i, rho_Arr(i));
    end
    mfprintf(f_a, '\n');
    
    mfprintf(f_a, 'Relaxation times\n');
    mfprintf(f_a, '  g: -\n');
    for i=1:length(ts_Arr)
        mfprintf(f_a, '  %d: %f\n', i, ts_Arr(i));
    end
    mfprintf(f_a, '\n');
    
    mclose(f_a);
end

//ЗАДАЁМ ФУНКЦИЮ ДЕЛЬТА f//WE SET THE DELTA FUNCTION f
//смотри между (44) и (45), (49)//look between (44) and (45), (49)
function y=delta_f(delta_f_wave, x, t, A)
    y=A*delta_f_wave*exp(%i*k_Const*x-w_Root*t);
    y=real(y);
endfunction

function X=uniform_Grid(a,b,N)//равномерная сетка на [a,b] с N узлами//uniform grid on [a,b] with N nodes
    h=(b-a)/(N-1);
    for n=1:N
        X(n)=a+h*(n-1);
    end 
endfunction

//РИСУЕМ ГРАФИКИ ПЛОТНОСТЕЙ И СКОРОСТЕЙ//DRAWING DENSITY AND VELOCITY GRAPHICS
if type_Of_Plot=='x' then
    X_x=uniform_Grid(0,С_Const,number_Of_Points);
    for i=1:length(X_x)
        Y_x_rho_g_t0(i)=rho_g+delta_f(delta_rho_g_Wave, X_x(i), 0, A_Density);
        Y_x_v_g_t0(i)=delta_f(delta_v_g_Wave, X_x(i), 0, A_Velocity);
        Y_x_rho_g(i)=rho_g+delta_f(delta_rho_g_Wave, X_x(i), t_Selected, A_Density);
        Y_x_v_g(i)=delta_f(delta_v_g_Wave, X_x(i), t_Selected, A_Velocity);
        
        if overwrite_Files=='on' then
            f_a=mopen(name_Arr_Double(1), 'at');
            mfprintf(f_a, '%d) %f         %f         %f\n', i, X_x(i), Y_x_rho_g_t0(i), Y_x_v_g_t0(i));
            mclose(f_a);
            
            f_a=mopen(name_Arr_Double(2+length(rho_Arr)), 'at');
            mfprintf(f_a, '%d) %f         %f         %f\n', i, X_x(i), Y_x_rho_g(i), Y_x_v_g(i));
            mclose(f_a);
        end
        
    end
    L=length(rho_Arr);
    
    subplot(L+1,2,1);
    plot2d(X_x, Y_x_rho_g_t0, style=color(color_Arr_Double(1)), strf='181');
    plot2d(X_x, Y_x_rho_g, style=color(color_Arr_Double(2)), strf='181');
    hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
    g=get('current_axes');
    g.title.text='Gas Density';
    
    subplot(L+1,2,2);
    plot2d(X_x, Y_x_v_g_t0, style=color(color_Arr_Double(1)), strf='181');
    plot2d(X_x, Y_x_v_g, style=color(color_Arr_Double(2)), strf='181');
    hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
    g=get('current_axes');
    g.title.text='Gas Velocity';
    
    for s=1:L
        for i=1:length(X_x)
            Y_x_rho_d_t0(i)=rho_Arr(s)+delta_f(delta_rho_d_Wave(s), X_x(i), 0, A_Density);
            Y_x_v_d_t0(i)=delta_f(delta_v_d_Wave(s), X_x(i), 0, A_Velocity);
            Y_x_rho_d(i)=rho_Arr(s)+delta_f(delta_rho_d_Wave(s), X_x(i), t_Selected, A_Density);
            Y_x_v_d(i)=delta_f(delta_v_d_Wave(s), X_x(i), t_Selected, A_Velocity);
            
            if overwrite_Files=='on' then
                f_a=mopen(name_Arr_Double(s+1), 'at');
                mfprintf(f_a, '%d) %f         %f         %f\n', i, X_x(i), Y_x_rho_d_t0(i), Y_x_v_d_t0(i));
                mclose(f_a);
                
                f_a=mopen(name_Arr_Double(s+2+L), 'at');
                mfprintf(f_a, '%d) %f         %f         %f\n', i, X_x(i), Y_x_rho_d(i), Y_x_v_d(i));
                mclose(f_a);
            end
            
        end
        
        subplot(L+1,2,2*s+1);
        plot2d(X_x, Y_x_rho_d_t0, style=color(color_Arr_Double(1)), strf='181');
        plot2d(X_x, Y_x_rho_d, style=color(color_Arr_Double(2)), strf='181');
        hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
        g=get('current_axes');
        g.title.text='Dust Density '+string(s);
        
        subplot(L+1,2,2*s+2);
        plot2d(X_x, Y_x_v_d_t0, style=color(color_Arr_Double(1)), strf='181');
        plot2d(X_x, Y_x_v_d, style=color(color_Arr_Double(2)), strf='181');
        hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
        g=get('current_axes');
        g.title.text='Dust Velocity '+string(s);
    end
elseif type_Of_Plot=='t' then
    X_t=uniform_Grid(0,T_Const,number_Of_Points);
    for i=1:length(X_t)
        Y_t_v_g(i)=delta_f(delta_v_g_Wave, x_Selected, X_t(i), 1)/Cs_Const;
        Y_t_rho_g(i)=delta_f(delta_rho_g_Wave, x_Selected, X_t(i), 1)/rho_g;
        
        if overwrite_Files=='on' then
            f_a=mopen(name_Arr_Normalized(1), 'at');
            mfprintf(f_a, '%d) %f       %f               %f\n', i, X_t(i), Y_t_rho_g(i), Y_t_v_g(i));
            mclose(f_a);
        end
        
    end
    
    legend_Arr(1)='Gas';
    
    subplot(2,1,1);
    plot2d(X_t, Y_t_v_g, style=color(color_Arr_Article(1)), strf='181');
    subplot(2,1,2);
    plot2d(X_t, Y_t_rho_g, style=color(color_Arr_Article(1)), strf='181');
    
    for i=1:length(rho_Arr)
        for j=1:length(X_t)
            Y_t_v_d(j)=delta_f(delta_v_d_Wave(i), x_Selected, X_t(j), 1)/Cs_Const;
            Y_t_rho_d(j)=delta_f(delta_rho_d_Wave(i), x_Selected, X_t(j), 1)/rho_Arr(i);
            
            if overwrite_Files=='on' then
                f_a=mopen(name_Arr_Normalized(i+1), 'at');
                mfprintf(f_a, '%d) %f       %f               %f\n', j, X_t(j), Y_t_rho_d(j), Y_t_v_d(j));
                mclose(f_a);
            end
            
        end
        
        legend_Arr(i+1)='Dust '+string(i);
        
        subplot(2,1,1);
        plot2d(X_t, Y_t_v_d, style=color(color_Arr_Article(i+1)), strf='181');
        hl=legend(legend_Arr, -1);
        g=get('current_axes');
        g.x_label.text='Time';
        g.y_label.text='Normalized Velocity';
        
        subplot(2,1,2);
        plot2d(X_t, Y_t_rho_d, style=color(color_Arr_Article(i+1)), strf='181');
        hl=legend(legend_Arr, -1);
        g=get('current_axes');
        g.x_label.text='Time';
        g.y_label.text='Normalized Density';
        
    end
end

//РАСЧЁТ ДЕЛЬТЫ ПЛОТНОСТЕЙ И СКОРОСТЕЙ//CALCULATIONS OF DENSITY AND VELOCITY DELTAS
//смотри (49)//see (49)
delta_rho_g=delta_f(delta_rho_g_Wave, x_Selected, 0, A_Density);
for i=1:length(delta_rho_d_Wave)
    delta_rho_i(i)=delta_f(delta_rho_d_Wave(i), x_Selected, 0, A_Density);
end
delta_v_g=delta_f(delta_v_g_Wave, x_Selected, 0, A_Velocity);
for i=1:length(delta_v_d_Wave)
    delta_v_i(i)=delta_f(delta_v_d_Wave(i), x_Selected, 0, A_Velocity);
end

//ЗАПИСЬ ДЕЛЬТЫ ПЛОТНОСТЕЙ И СКОРОСТЕЙ В ФАЙЛ//DENSITY AND VELOCITY DELTAS RECORDING TO FILE
if overwrite_Files=='on' then
    f_a=mopen(full_Name, 'at');
    mfprintf(f_a, 'Density deltas at t=0\n');
    sAD=string(A_Density);
    CN=delta_rho_g_Wave;
    mfprintf(f_a, '  g: %s*(%f*cos(%f*x) %s %f*sin(%f*x))\n',sAD, real(CN), k_Const, imag_Sign_Rev(CN), abs(imag(CN)), k_Const);
    for i=1:length(delta_rho_i)
        CN=delta_rho_d_Wave(i);
        mfprintf(f_a, '  %d: %s*(%f*cos(%f*x) %s %f*sin(%f*x))\n', i, sAD, real(CN), k_Const, imag_Sign_Rev(CN), abs(imag(CN)), k_Const);
    end
    mfprintf(f_a, '\n');
    
    mfprintf(f_a, 'Velocity deltas at t=0\n');
    sAV=string(A_Velocity);
    CN=delta_v_g_Wave;
    mfprintf(f_a, '  g: %s*(%f*cos(%f*x) %s %f*sin(%f*x))\n',sAV, real(CN), k_Const, imag_Sign_Rev(CN), abs(imag(CN)), k_Const);
    for i=1:length(delta_v_i)
        CN=delta_v_d_Wave(i);
        mfprintf(f_a, '  %d: %s*(%f*cos(%f*x) %s %f*sin(%f*x))\n', i, sAV, real(CN), k_Const, imag_Sign_Rev(CN), abs(imag(CN)), k_Const);
    end
    mclose(f_a);
end
