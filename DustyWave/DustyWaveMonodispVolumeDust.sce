//НАЧАЛЬНЫЕ ДАННЫЕ
t_stop=10000000;//время релаксации
rho_g=1;//плотность газа
rho_d=0.01;//плотность пыли
k=1;//количество волн
a_g=1;//скорость звука в газе
a_d=0;//скорость звука в пыли
theta=0.4;//объёмная доля

delta_rho_g_wave=0.01;

С_Const=1.0;//длина отрезка координат
t_Selected=0.02;//в какой момент времени рассматриваем волну
number_Of_Points=100;//количество узлов, по которым строим графики

overwrite_Files='on';//перезапись файлов при отрисовке графиков (on-да)

path='D:\Projects\Multigrain\DispersionRel\ScilabSovlers';//путь до текстовых файлов
//path='C:\Users\User\Desktop\ГАЗОПЫЛЬ\program_for_article\Output files';//путь до текстовых файлов

//КОНЕЦ НАЧАЛЬНЫХ ДАННЫХ

k=(2*%pi*k)/С_Const;
epsilon=rho_d/rho_g;

//КОЭФФИЦИЕНТЫ ПОЛИНОМА
S(1)=k^4*(a_g^2*(a_d^2+theta^2*a_g^2*rho_g/(rho_d*(1-theta)^2))-a_g^2*a_d^2*theta^2/(1-theta)^2);
S(2)=-k^2/t_stop*(theta*(a_d^2+a_g^2)/(1-theta)+a_g^2+rho_d/rho_g*a_d^2+a_g^2*theta^2/(1-theta)^2);
S(3)=k^2*(a_d^2+(rho_g*theta^2*a_g^2)/(rho_d*(1-theta)^2)+a_g^2);
S(4)=-(1+rho_d/rho_g)/t_stop;
S(5)=1;

//КОРНИ ПОЛИНОМА
for i=1:length(S)
    coef(i)=S(6-i);
end
R=roots(coef);
disp(R);

//КОРЕНЬ С ОТРИЦАТЕЛЬНОЙ МНИМОЙ ЧАСТЬЮ
omega=0;
for i=1:length(R)
    if imag(R(i))>6.3 then
        omega=R(i);
    end
end
disp(omega);

//ФУНКЦИИ ДЛЯ РАСЧЁТА ДЕЛЬТ
function y=delta_rho_g(x,t)
    y=delta_rho_g_wave*exp(%i*k*x-omega*t);
endfunction

function y=delta_u_g(x,t)
    y=delta_rho_g_wave*exp(%i*k*x-omega*t);
    y=-%i*omega/(k*rho_g)*y;
endfunction

function y=delta_rho_d(x,t)
    y=delta_rho_g_wave*exp(%i*k*x-omega*t);
    y=((-omega+epsilon/t_stop)*omega/(k*rho_g)-k*a_g^2/rho_g)*y;
    y=y/(k*a_g^2/(rho_d*(1/theta-1))+epsilon*omega/(t_stop*k*rho_d));
endfunction

function y=delta_u_d(x,t)
    y=delta_rho_g_wave*exp(%i*k*x-omega*t);
    y=((-omega+epsilon/t_stop)*omega/(k*rho_g)-k*a_g^2/rho_g)*y;
    y=y/(k*a_g^2/(rho_d*(1/theta-1))+epsilon*omega/(t_stop*k*rho_d));
    y=-%i*omega/(k*rho_d)*y;
endfunction

//РАСЧЁТ КОЭФФИЦИЕНТОВ ПЕРЕД КОСИНУСОМ В ДЕЛЬТАХ
coef_cos_delta_rho_g=real(delta_rho_g_wave);

coef_cos_delta_u_g=real(-%i*omega/(k*rho_g)*delta_rho_g_wave);

coef_cos_delta_rho_d=real(((-omega+epsilon/t_stop)*omega/(k*rho_g)-k*a_g^2/rho_g)*delta_rho_g_wave/(k*a_g^2/(rho_d*(1/theta-1))+epsilon*omega/(t_stop*k*rho_d)));

coef_cos_delta_u_d=real(-%i*omega/(k*rho_d)*((-omega+epsilon/t_stop)*omega/(k*rho_g)-k*a_g^2/rho_g)*delta_rho_g_wave/(k*a_g^2/(rho_d*(1/theta-1))+epsilon*omega/(t_stop*k*rho_d)));

//РАСЧЁТ КОЭФФИЦИЕНТОВ ПЕРЕД СИНУСОМ В ДЕЛЬТАХ
coef_sin_delta_rho_g=imag(delta_rho_g_wave);

coef_sin_delta_u_g=imag(-%i*omega/(k*rho_g)*delta_rho_g_wave);

coef_sin_delta_rho_d=imag(((-omega+epsilon/t_stop)*omega/(k*rho_g)-k*a_g^2/rho_g)*delta_rho_g_wave/(k*a_g^2/(rho_d*(1/theta-1))+epsilon*omega/(t_stop*k*rho_d)));

coef_sin_delta_u_d=imag(-%i*omega/(k*rho_d)*((-omega+epsilon/t_stop)*omega/(k*rho_g)-k*a_g^2/rho_g)*delta_rho_g_wave/(k*a_g^2/(rho_d*(1/theta-1))+epsilon*omega/(t_stop*k*rho_d)));

//РЕАЛЬНАЯ ЧАСТЬ ДЕЛЬТЫ В ОПРЕДЕЛЁННЫЙ МОМЕНТ ВРЕМЕНИ
function y=real_delta(delta,x,t)
    y=real(delta(x,t));
endfunction

//ФУНКЦИЯ ДЛЯ ПОСТРОЕНИЯ СЕТКИ
function X=uniform_Grid(a,b,N)//равномерная сетка на [a,b] с N узлами
    h=(b-a)/(N-1);
    for n=1:N
        X(n)=a+h*(n-1);
    end 
endfunction

X=uniform_Grid(0,С_Const,number_Of_Points);
for i=1:length(X)
    RDRG0(i)=real_delta(delta_rho_g,X(i),0);
    RDUG0(i)=real_delta(delta_u_g,X(i),0);
    RDRD0(i)=real_delta(delta_rho_d,X(i),0);
    RDUD0(i)=real_delta(delta_u_d,X(i),0);
    
    RDRGt(i)=real_delta(delta_rho_g,X(i),t_Selected);
    RDUGt(i)=real_delta(delta_u_g,X(i),t_Selected);
    RDRDt(i)=real_delta(delta_rho_d,X(i),t_Selected);
    RDUDt(i)=real_delta(delta_u_d,X(i),t_Selected);
end

//СОЗДАЕМ ПУСТОЙ ТЕКСТОВЫЙ ФАЙЛ
if overwrite_Files=='on' then
    full_Name=path+'\t=0.txt';
    f_w=mopen(full_Name, 'wt');
    mfprintf(f_w, '#real(delta_rho_g)(x)=('+string(coef_cos_delta_rho_g)+')*cos('+string(k)+'*x)-('+string(coef_sin_delta_rho_g)+')*sin('+string(k)+'*x)\n');
    mfprintf(f_w, '#real(delta_u_g)(x)=('+string(coef_cos_delta_u_g)+')*cos('+string(k)+'*x)-('+string(coef_sin_delta_u_g)+')*sin('+string(k)+'*x)\n');
    mfprintf(f_w, '#real(delta_rho_d)(x)=('+string(coef_cos_delta_rho_d)+')*cos('+string(k)+'*x)-('+string(coef_sin_delta_rho_d)+')*sin('+string(k)+'*x)\n');
    mfprintf(f_w, '#real(delta_u_d)(x)=('+string(coef_cos_delta_u_d)+')*cos('+string(k)+'*x)-('+string(coef_sin_delta_u_d)+')*sin('+string(k)+'*x)\n');
    mfprintf(f_w, '\n');
    mfprintf(f_w, '#      x          real(delta_rho_g)  real(delta_u_g)  real(delta_rho_d)  real(delta_u_d)\n');
    for i=1:length(X)
        mfprintf(f_w, '%d) %f          %f          %f          %f          %f\n', i, X(i), RDRG0(i), RDUG0(i), RDRD0(i), RDUD0(i));
    end
    mclose(f_w);
end

if overwrite_Files=='on' then
    full_Name=path+'\t='+string(t_Selected)+'.txt';
    f_w=mopen(full_Name, 'wt');
    mfprintf(f_w, '#      x          real(delta_rho_g)  real(delta_u_g)  real(delta_rho_d)  real(delta_u_d)\n');
    for i=1:length(X)
        mfprintf(f_w, '%d) %f          %f          %f          %f          %f\n', i, X(i), RDRGt(i), RDUGt(i), RDRDt(i), RDUDt(i));
    end
    mclose(f_w);
end

//СТРОИМ ГРАФИКИ
subplot(4,1,1);
plot(X,RDRG0,'blue');
plot(X,RDRGt,'green');
hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
g=get('current_axes');
g.title.text='real(delta_rho_g)';
    
subplot(4,1,2);
plot(X,RDUG0,'blue');
plot(X,RDUGt,'green');
hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
g=get('current_axes');
g.title.text='real(delta_u_g)';

subplot(4,1,3);
plot(X,RDRD0,'blue');
plot(X,RDRDt,'green');
hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
g=get('current_axes');
g.title.text='real(delta_rho_d)';

subplot(4,1,4);
plot(X,RDUD0,'blue');
plot(X,RDUDt,'green');
hl=legend(['t = 0', 't = '+string(t_Selected)], -1);
g=get('current_axes');
g.title.text='real(delta_u_d)';
