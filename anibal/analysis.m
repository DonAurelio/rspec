% Accelerograma
function obj = Acelerograma(signal, sensivility, fact,fm) 
    obj.sampleRate=fm;            
    obj.volts  = detrend((signal));            
%   obj.gal(:,1)=detrend(obj.volts(:,1));
%   obj.gal(:,2)=detrend(obj.volts(:,2));
%   obj.gal(:,3)=detrend(obj.volts(:,3));
    obj.T=1/fm;
    obj.timeTotal=(0:obj.T:(obj.T*length(signal))-obj.T)';
    obj.g   = (obj.volts./1.25);
%   aux=obj.g.*981;
    obj.gal = obj.g.*981;  
    
end


% RAC_GUI
%% espectro de respuesta
% variables
ro=0.05; % Amortiguacion
T_pas=0:0.01:4; % Paso para el periodo
wn_pas=2*pi./T_pas; % Frecuencia

[station] = espectroRespuesta(station,ro,wn_pas,T_pas);


% Método que permite calcular el espectro de respuesta
function [obj] = espectroRespuesta(obj,z,wn,Tp)
           
     RsueE=zeros(1,length(wn)); 
     RsueE2=zeros(1,length(wn));
     RsueE3=zeros(1,length(wn));
     for i=1:length(wn)
         H=tf(-1,[1 2*z*wn(i) wn(i).^2]);

         acesE=lsim(H,(obj.gal(:,1)),obj.timeTotal);
         acesE2=lsim(H,(obj.gal(:,2)),obj.timeTotal);
         acesE3=lsim(H,(obj.gal(:,3)),obj.timeTotal);
        % Respuesta maxima
     RsueX(i)=max(acesE); 
     RsueY(i)=max(acesE2);  
     RsueZ(i)=max(acesE3);
     end
    obj.esp_respuesta(:,1)=RsueX'.*(wn.^2)';
    obj.esp_respuesta(:,2)=RsueY'.*(wn.^2)';
    obj.esp_respuesta(:,3)=RsueZ'.*(wn.^2)';
    
    respuesta=figure('Name','Espectro Respuesta','Position',[10 10 650 250])
    subplot(1,3,1)
    plot(Tp,obj.esp_respuesta(:,1))
    ylabel ('Aceleración (cm/seg/seg)')
    axis square;
    title ('N-S'); 

    subplot(1,3,2)
    plot(Tp,obj.esp_respuesta(:,2))
    xlabel ('Tiempo (seg)')
    axis square;
    title ('E-O'); 

    subplot(1,3,3)
    plot(Tp,obj.esp_respuesta(:,3))
    axis square;
    title ('Ver');
    
%             F = getframe(respuesta);
%             imwrite(F.cdata,'Respuesta.jpg');
end 