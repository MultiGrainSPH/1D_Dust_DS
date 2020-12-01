//НАЧАЛЬНЫЕ ДАННЫЕ//INITIAL DATA
rho_g=1;//плотность газа//gas density
rho_d=1;//плотность пыли//dust density
t_stop=0.000001;//время релаксации//relaxation time
k=1;//количество волн//number of waves
c_g=1;//скорость звука в газе//speed of sound in gas
c_d=0.0;//скорость звука в пыли//speed of sound in dust
theta=0.0;//объёмная доля//volume fraction

delta_rho_g_wave=0.01;

С_Const=1;//длина отрезка координат//coordinate length
t_Selected=0.3;//в какой момент времени рассматриваем волну//at what point in time we consider the wave
number_Of_Points=100;//количество узлов, по которым строим графики//the number of nodes on which to build graphs

overwrite_Files='on';//перезапись файлов при отрисовке графиков (on-да)//overwrite files when drawing graphs (on-yes)

path='D:\Projects\Multigrain\DispersionRel\ScilabSovlers';//путь до текстовых файлов//path to text files
//path='C:\Users\User\Desktop\ГАЗОПЫЛЬ\MonoDisp\Output files';//путь до текстовых файлов

//КОНЕЦ НАЧАЛЬНЫХ ДАННЫХ//END OF INITIAL DATA

k=(2*%pi*k)/С_Const;
epsilon=rho_d/rho_g;
xi=theta/(1-theta);

//КОЭФФИЦИЕНТЫ ПОЛИНОМА//POLYNOMIC COEFFICIENTS
S(1)=k^4*(c_g^2*c_d^2);
S(2)=-k^2/t_stop*((xi+1)^2*c_g^2+epsilon*c_d^2);
S(3)=k^2*((xi^2/epsilon+1)*c_g^2+c_d^2);
S(4)=-(1+epsilon)/t_stop;
S(5)=1;

//КОРНИ ПОЛИНОМА//ROOTS OF POLINOMA
for i=1:length(S)
    coef(i)=S(6-i);
end
R=roots(coef);
disp(R);

//КОРЕНЬ С ПОЛОЖИТЕЛЬНОЙ МНИМОЙ ЧАСТЬЮ//ROOT WITH A POSITIVE IMAGINARY PART
omega=0;
for i=1:length(R)
    if imag(R(i))>0 then
        omega=R(i);
    end
end
disp(omega);

//ДЕЛЬТЫ С ВОЛНАМИ//
delta_u_g_wave = -%i*omega/(k*rho_g)*delta_rho_g_wave;
delta_rho_d_wave = ((-omega+epsilon/t_stop)*omega/(k*rho_g)-k*c_g^2/rho_g)/(k*c_g^2/(rho_d*(1/theta-1))+epsilon*omega/(t_stop*k*rho_d))*delta_rho_g_wave;
delta_u_d_wave = -%i*omega/(k*rho_g)*delta_rho_d_wave;

//ФУНКЦИИ ДЛЯ РАСЧЁТА ДЕЛЬТ//FUNCTIONS FOR CALCULATING DELTA
function y=delta_rho_g(x,t)
    y=delta_rho_g_wave*exp(%i*k*x-omega*t);
endfunction

function y=delta_u_g(x,t)
    y=delta_u_g_wave*exp(%i*k*x-omega*t);
endfunction

function y=delta_rho_d(x,t)
    y=delta_rho_d_wave*exp(%i*k*x-omega*t);
endfunction

function y=delta_u_d(x,t)
    y=delta_u_d_wave*exp(%i*k*x-omega*t);
endfunction

//РАСЧЁТ КОЭФФИЦИЕНТОВ ПЕРЕД КОСИНУСОМ В ДЕЛЬТАХ//CALCULATION OF COEFFICIENTS BEFORE THE COSINUS IN DELTA
coef_cos_delta_rho_g=real(delta_rho_g_wave);

coef_cos_delta_u_g=real(delta_u_g_wave);

coef_cos_delta_rho_d=real(delta_rho_d_wave);

coef_cos_delta_u_d=real(delta_u_d_wave);

//РАСЧЁТ КОЭФФИЦИЕНТОВ ПЕРЕД СИНУСОМ В ДЕЛЬТАХ//CALCULATION OF COEFFICIENTS BEFORE SINUS IN DELTA
coef_sin_delta_rho_g=imag(delta_rho_g_wave);

coef_sin_delta_u_g=imag(delta_u_g_wave);

coef_sin_delta_rho_d=imag(delta_rho_d_wave);

coef_sin_delta_u_d=imag(delta_u_d_wave);

//ФУНКЦИЯ ДЛЯ ПОСТРОЕНИЯ СЕТКИ//FUNCTION FOR CONSTRUCTION GRID
function X=uniform_Grid(a,b,N)//равномерная сетка на [a,b] с N узлами//uniform grid on [a, b] with N nodes
    h=(b-a)/(N-1);
    for n=1:N
        X(n)=a+h*(n-1);
    end 
endfunction

X=uniform_Grid(0,С_Const,number_Of_Points);
for i=1:length(X)
    RDRG0(i)=real(delta_rho_g(X(i),0));
    RDUG0(i)=real(delta_u_g(X(i),0));
    RDRD0(i)=real(delta_rho_d(X(i),0));
    RDUD0(i)=real(delta_u_d(X(i),0));
    
    RDRGt(i)=real(delta_rho_g(X(i),t_Selected));
    RDUGt(i)=real(delta_u_g(X(i),t_Selected));
    RDRDt(i)=real(delta_rho_d(X(i),t_Selected));
    RDUDt(i)=real(delta_u_d(X(i),t_Selected));
end

//СОЗДАЕМ ПУСТОЙ ТЕКСТОВЫЙ ФАЙЛ//CREATING AN EMPTY TEXT FILE
if overwrite_Files=='on' then
    full_Name=path+'\t=0.txt';
    f_w=mopen(full_Name, 'wt');
    mfprintf(f_w, '#real_delta_rho_g(x)=('+string(coef_cos_delta_rho_g)+')*cos('+string(k)+'*x)-('+string(coef_sin_delta_rho_g)+')*sin('+string(k)+'*x)\n');
    mfprintf(f_w, '#real_delta_u_g(x)=('+string(coef_cos_delta_u_g)+')*cos('+string(k)+'*x)-('+string(coef_sin_delta_u_g)+')*sin('+string(k)+'*x)\n');
    mfprintf(f_w, '#real_delta_rho_d(x)=('+string(coef_cos_delta_rho_d)+')*cos('+string(k)+'*x)-('+string(coef_sin_delta_rho_d)+')*sin('+string(k)+'*x)\n');
    mfprintf(f_w, '#real_delta_u_d(x)=('+string(coef_cos_delta_u_d)+')*cos('+string(k)+'*x)-('+string(coef_sin_delta_u_d)+')*sin('+string(k)+'*x)\n');
    mfprintf(f_w, '\n');
    mfprintf(f_w, '#      x          real(delta_rho_g)  real(delta_u_g)  real(delta_rho_d)  real(delta_u_d)\n');
    for i=1:length(X)
        mfprintf(f_w, '%f          %f          %f          %f          %f\n', X(i), RDRG0(i), RDUG0(i), RDRD0(i), RDUD0(i));
    end
    mclose(f_w);
end

if overwrite_Files=='on' then
    full_Name=path+'\t='+string(t_Selected)+'.txt';
    f_w=mopen(full_Name, 'wt');
    mfprintf(f_w, '#      x          real(delta_rho_g)  real(delta_u_g)  real(delta_rho_d)  real(delta_u_d)\n');
    for i=1:length(X)
        mfprintf(f_w, '%f          %f          %f          %f          %f\n', X(i), RDRGt(i), RDUGt(i), RDRDt(i), RDUDt(i));
    end
    mclose(f_w);
end

//СТРОИМ ГРАФИКИ//BUILDING GRAPHICS
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
