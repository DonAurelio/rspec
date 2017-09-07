%% 
clc, clear, close all

y=load('panamer.txt');
y=y(1:length(y));
y1=load('Resultado_Capa1.txt');
y2=load('Seudo_Reforma.txt');
t2=y2(: ,1);
y2=y2(: , 2);
% y1=load('Resultado2_Capa1.txt');
y1=y1(1:length(y1))/9.81;
Ts=.005;
% Ts1=0.02;
Ts1=0.005;

%y= y/(max(abs(y))/0.25);
Fs=1/Ts;
Fs1=1/Ts1;
t=(0:length(y)-1)/Fs;
t1=(0:length(y1)-1)/Fs1;
% [N D] = butter(2,5/Fs/2, 'high');
% yfilt= filtfilt(N,D,y);

subplot(3,2,1); %Plot Panamericano
plot(t,y,'r'); grid on,hold on
%  axis([0 max(t)  -0.3 0.3]);
 ylabel('Entrada en aceleración_P [g]');
 xlabel('Tiempo [s]');

subplot(3,2,2); % Plot Resultado
plot(t1,y1,'r'); grid on,hold on
%  axis([0 max(t1)  -0.3 0.3]);
 ylabel('Entrada en aceleración_R [g]');
 xlabel('Tiempo [s]');

%Paso a dominio de frecuencia  Panamericana y plot

Y=abs(fft(y));
f=linspace(0,Fs,length(y));
subplot (3,2,3)
plot (f,2*Y/length(y)); 
xlim([0 Fs/2]) % Hice cambio para observar mejor la grafica //ademas del filtro para eliminar la primera frecuencia
xlabel('Frecuencia_P [s]');

%Paso a dominio de frecuencia  Respuesta y plot
Y1=abs(fft(y1));
f1=linspace(0,Fs1,length(y1));
subplot (3,2,4)
plot (f1,2*Y1/length(y1)); 
xlim([0 Fs1/2]) % Hice cambio para observar mejor la grafica //ademas del filtro para eliminar la primera frecuencia
xlabel('Frecuencia_R [s]');

%Espectro de respuesta en seudo aceleracion Panamericana
subplot(3,2,5)
Tn=.001:.01:10;
PSA=funcesp(y,0.05,Fs,Tn)';
plot(Tn,PSA);
 ylabel('Pseudo_Aceleración_P [g]');
 xlabel('Periodo_P [s]');


%Espectro de respuesta en seudo aceleracion resultado
subplot(3,2,6)
Tn1=.001:.01:10;
PSA1=funcesp(y1,0.05,Fs1,Tn1)';
plot(Tn1,PSA1), hold on
plot(t2,y2);
 ylabel('Pseudo_Aceleración_R [g]');
 xlabel('Periodo_R [s]');

 fita1=(100*(1-(norm(abs()-abs(ED')))/(norm(abs(EP)-mean(EP)))));
 %%
% 
% 
% for n=1: length(Tn)
% T=Tn(n);
% Z=0.05;
% wn=2*pi/T;
% 
% H=tf([0 0 -1],[1 2*Z*wn wn^2]);
% H=c2d(H,1/Fs);
% D= lsim(H,y*9.81,t); %% se multiplica por 9.81 por que se debe pasar a unidade de desplazamiento
% H2=tf([-1 0 0],[1 2*Z*wn wn^2]);%% se multiplica por 9.81 por que se debe pasar a unidades a aceleraccion
% H2=c2d(H2,1/Fs);
% A= lsim(H2,y*9.81,t)/9.81;
% 
% % subplot(3,2,3);
% % plot(t,D); grid on
% % ylabel('Resp. en desplazamiento [m]');
% % xlabel('Tiempo [s]');
% 
% % subplot(3,2,5);
% % plot(t,A); grid on
% % ylabel('Resp. en aceleración [g]');
% % xlabel('Tiempo [s]');
% 
% P(n)= max(abs(D));
% 
% % subplot(2,2,3);
% % plot(Tn(1:n),P),grid on
% % ylabel ('Esp. respuesta desp. [m]');
% % xlabel ('Periodo [s]');
% 
% SA(n)=P(n)*(2*pi/T)^2/9.81;
% % subplot(2,2,4);
% % plot(Tn(1:n),SA),grid on, hold on
% % ylabel ('Esp. respuesta Seudo- aceleracion [g]');
% % xlabel ('Periodo [s]');
% 
% A1(n)=max(abs(A))/9.81;
% % subplot(2,2,4);
% % plot(Tn(1:n),A1),grid on
% % legend('Esp. respuesta Seudo-aceleración [g]','Esp. respuesta aceleración [g]')
% %  ylabel ('Esp. respuesta Aceleracion [g]');
% %  xlabel ('Periodo [s]');
% 
% end
%%