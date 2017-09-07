%% ESPECTRO DE RESPUESTA
tic
clc, clear all, warning off
t=0.02:0.02:4; % Vector de periodos para graficar
T=1;
z=[0.05 0.1 0.2]; % Razon de amortiguamiento
m=1; % Masa
fm=200;% Frecuencia de muestreo del sensor
sv=10; % Sensibilidad en voltios/gravedad
gal=981;% Gravedad en gales
fe=1*gal/sv;% Factor de escala
d=1/fm; % Delta tiempo
reg=load('prueba1.txt');  
% ti= 0:d:(length(reg)-1)*d;%  reg(:,1);%0:0.1:1;    
    fs=detrend(reg(1:8000,1),'constant') *fe;
    fx=[0.1 50];
    den=((1/d)/2);
    vf=[fx(1)/den fx(2)/den];
    [C1,C2]=butter(2,vf,'bandpass');
    fst=filtfilt(C1,C2,fs);
    y=  fst ;
for g=1:3
for i=1:3
     
    for j=1:length(t)
       
    Wn=(2*pi)/t(j);
    k=m*(Wn^2);
    a1=(2.71828)^(-z(g)*Wn*d);
    Wd=Wn*sqrt(1-(z(g)^2));
    a2=sin(Wd*d);
    a3=cos(Wd*d);
    
    A=a1*(((z(g)/(sqrt(1-z(g)^2)))*a2)+a3);
    B=a1*((1/Wd)*a2);
    C=(1/k)*((2*z(g)/(Wn*d))+a1*((((1-2*(z(g)^2))/(Wd*d))-(z(g)/(sqrt(1-z(g)^2))))*a2-(1+(2*z(g)/(Wn*d)))*a3));
    D=(1/k)*(1-(2*z(g)/(Wn*d))+a1*(((2*(z(g)^2)-1)/(Wd*d))*a2+(2*z(g)/(Wn*d))*a3));
    A1=-a1*(((Wn/(sqrt(1-z(g)^2)))*a2));
    B1=a1*(a3-((z(g)/(sqrt(1-z(g)^2)))*a2));
    C1=(1/k)*((-1/d)+a1*(((Wn/(sqrt(1-(z(g)^2))))+(z(g)/(d*sqrt(1-(z(g)^2)))))*a2+(a3/d)));
    D1=(1/(k*d))*(1-a1*((z(g)/(sqrt(1-(z(g)^2))))*a2+a3));
   
         
    ui=[];
    upi=[];   
    ui(1)=0;
    upi(1)=0;
    
    for k=1:(length(y)-1);
    ui(k+1)=A*ui(k)+B*upi(k)+C*y(k)+D*y(k+1);
    upi(k+1)=A1*ui(k)+B1*upi(k)+C1*y(k)+D1*y(k+1);
    
    end
    
    Dmax(j,i,g)=max(abs(ui));
    Vmax(j,i,g)=max(abs(ui*Wn));
    Amax(j,i,g)=max(abs((ui*(Wn^2))));
    
            
    end
    
end

DA(1:length(t),g)=mean(Amax(:,:,g),2);
DV(1:length(t),g)=mean(Vmax(:,:,g),2);
DD(1:length(t),g)=mean(Dmax(:,:,g),2);
end


close all
for i=1

    figure (i);
    subplot(3,1,1)
     hold on
    plot(t,Dmax(:,i,1),t,Dmax(:,i,2),t,Dmax(:,i,3))
    semilogx(t,Dmax(:,i,g))
    legend('z=5%','z=10%','z=20%');
    xlabel('Tn [s]'),ylabel('D [m]')
    %xlim([0,10])
    title('Espectro de Respuesta en Desplazamiento')
    subplot(3,1,2)
     hold on
    plot(t,Vmax(:,i,1),t,Vmax(:,i,2),t,Vmax(:,i,3))
    semilogx(t,Vmax(:,i,g))
    legend('z=5%','z=10%','z=20%');
    xlabel('Tn [s]'),ylabel('V [m/s]')
    %xlim([0,10])
    title('Espectro de Respuesta en Pseudovelocidad')
    subplot(3,1,3)
     hold on
    plot(t,Amax(:,i,1),t,Amax(:,i,2),t,Amax(:,i,3))
    legend('z=5%','z=10%','z=20%');
    semilogx(t,Amax(:,i,g))
    xlabel('Tn [s]'),ylabel('A [cm/s^2]')
    %xlim([0,10])
    title('Espectro de Respuesta en Pseudoaceleración')   
    
end


%%
f2=zeros(3000,3);
f4=zeros(3000,3);
f6=zeros(3000,3);
for g=1:3
   for i=1
       for h=1:3000
       f1(h,i,g)=((Amax(h,i,g)-DA(h,g))^2)/10;
       f2(h,g)=f2(h,g)+f1(h,i,g);
       SA(h,g)=sqrt(f2(h,g));
       f3(h,i,g)=((Dmax(h,i,g)-DD(h,g))^2)/10;
       f4(h,g)=f4(h,g)+f3(h,i,g);
       SD(h,g)=sqrt(f4(h,g));
       f5(h,i,g)=((Vmax(h,i,g)-DV(h,g))^2)/10;
       f6(h,g)=f6(h,g)+f5(h,i,g);
       SV(h,g)=sqrt(f6(h,g));
       end
       
   end   
    
    
end
%%
for i=1:3000
    if t(i)<=0.7
        sa(i)=2.5*0.25*1.3;
    end
    if t(i)>0.7 && t(i)<=4.56
        sa(i)=1.2*0.25*1.9/t(i);
    end
    if t(i)>4.56
        sa(i)=(1.2*0.25*1.9*4.56)/t(i)^2;
    end
end
figure(45);%'Name',['Espectros de Diseño con zita=',num2str(z*100)'','%'],'NumberTitle','off','Color',[0.6 0.8 0.89])

% subplot(3,1,1)
% plot(t,DD(:,1)+SD(:,1))
% legend('z=5%')%,'z=10%','z=20%');
% %loglog(t,DD)
% xlabel('Tn [s]'),ylabel('D[m]')
% xlim([0,5])
% title('Espectro de Diseño en Desplazamiento')
% subplot(3,1,2)
% plot(t,DV+SV)
% legend('z=5%','z=10%','z=20%');
% %loglog(t,DV)
% xlabel('Tn [s]'),ylabel('V [m/s]')
% xlim([0,5])
% title('Espectro de Diseño en Pseudovelocidad')
% subplot(3,1,3)
plot(t,DA(:,1)+SA(:,1),t,sa)
legend('Real','Calculado según la norma')%,'z=10%','z=20%');
%loglog(t,DA)
xlabel('Tn [s]'),ylabel('A [g]')
xlim([0,10])
title('Espectro de Diseño en Pseudoaceleracion para z=5%')


%%
%desviación estandar
f2=zeros(3000,3);
for g=1
   for i=1:10
       for h=1:3000
       f1(h,i,g)=((Amax(h,i,g)-DA(h,g))^2)/10;
       f2(h,g)=f2(h,g)+f1(h,i,g);
       S(h,g)=sqrt(f2(h,g));
       end
       
   end   
    
    
end
figure(35)
plot(t,S)

%%
clc, clear all

nn=11; % nodos en planta
p=5;%pisos
ngdl=nn*(p+1)*6;
mr=5;%# de marcos
ne=(2*nn+mr-1)*p;

