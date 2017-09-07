%% ESPECTRO DE RESPUESTA
tic
clc, clear all, warning off
%pkg install -forge control signal 
pkg load signal

%% Reading comand line parameters
args = argv();
file = strvcat(args(1,1));
damping = str2num(strvcat(args(2,1)));
sample_rate = str2num(strvcat(args(3,1)));
reg = load(file);

disp(file)
disp(damping)
disp(sample_rate)

t=0.02:0.02:4; % Vector de periodos de oscilacion
T=1;
z=damping; % Razon de amortiguamiento (0.05)
m=1; % Masa
fm=sample_rate;% Frecuencia de muestreo del sensor (200 MHz)
sv=10; % Sensibilidad en voltios/gravedad
gal=981;% Gravedad en gales
fe=1*gal/sv;% Factor de escala
d=1/fm; % Delta tiempo

fs=detrend(reg(1:length(reg),1),'constant') *fe;
fx=[0.1 50];
den=((1/d)/2);
vf=[fx(1)/den fx(2)/den];
[C1,C2]=butter(2,vf,'pass');
fst=filtfilt(C1,C2,fs);
y=fst;

for j=1:length(t)
   
    Wn=(2*pi)/t(j);
    k=m*(Wn^2);
    a1=(2.71828)^(-z*Wn*d);
    Wd=Wn*sqrt(1-(z^2));
    a2=sin(Wd*d);
    a3=cos(Wd*d);

    A=a1*(((z/(sqrt(1-z^2)))*a2)+a3);
    B=a1*((1/Wd)*a2);
    C=(1/k)*((2*z/(Wn*d))+a1*((((1-2*(z^2))/(Wd*d))-(z/(sqrt(1-z^2))))*a2-(1+(2*z/(Wn*d)))*a3));
    D=(1/k)*(1-(2*z/(Wn*d))+a1*(((2*(z^2)-1)/(Wd*d))*a2+(2*z/(Wn*d))*a3));
    A1=-a1*(((Wn/(sqrt(1-z^2)))*a2));
    B1=a1*(a3-((z/(sqrt(1-z^2)))*a2));
    C1=(1/k)*((-1/d)+a1*(((Wn/(sqrt(1-(z^2))))+(z/(d*sqrt(1-(z^2)))))*a2+(a3/d)));
    D1=(1/(k*d))*(1-a1*((z/(sqrt(1-(z^2))))*a2+a3));

         
    ui=[];
    upi=[];   
    ui(1)=0;
    upi(1)=0;

    for k=1:(length(y)-1);
        ui(k+1)=A*ui(k)+B*upi(k)+C*y(k)+D*y(k+1);
        upi(k+1)=A1*ui(k)+B1*upi(k)+C1*y(k)+D1*y(k+1);
    endfor

Dmax(j,1)=max(abs(ui));
Vmax(j,1)=max(abs(ui*Wn));
Amax(j,1)=max(abs((ui*(Wn^2))));

endfor

save ("-ascii","dmax.out", "Dmax");
save ("-ascii","vmax.out", "Vmax");
save ("-ascii","amax.out", "Amax");