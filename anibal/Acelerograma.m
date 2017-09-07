classdef Acelerograma

    %% Propiedades generales:
    properties
        % Frecuencia de muestreo
        sampleRate
        
        % señal en voltios
        volts 
        
        % señal en g
        g

        % señal en gales
        gal
        
        % tiempo total de la señal
        timeTotal
        
        % periodo de la señal
        T
        
        %espectro de fourier
        fourier
        
        %frecuencias de fourier
        fourier_HZ
        
        %espectro respuesta
        esp_respuesta
        
        %maximos
        maxs
        
    end
    
    %% Métodos que requieren argumentos
    methods
        %% Constructor
        function obj = Acelerograma(signal, sensivility, fact,fm) 
            obj.sampleRate=fm;            
            obj.volts  = detrend((signal));            
%             obj.gal(:,1)=detrend(obj.volts(:,1));
%             obj.gal(:,2)=detrend(obj.volts(:,2));
%             obj.gal(:,3)=detrend(obj.volts(:,3));
            obj.T=1/fm;
            obj.timeTotal=(0:obj.T:(obj.T*length(signal))-obj.T)';
            obj.g   = (obj.volts./1.25);
%             aux=obj.g.*981;
            obj.gal = obj.g.*981;  
            
        end

        %% Metodos        
        %% periodo
        % Permite actualizar el periodo
        function obj = setT(obj, T)
            obj.ro = T;
        end
      
        % Obtiene el periodo
        function get_ro = getT(obj)
            get_ro = obj.T;
        end
        
        %% señal en Volts
        % Permite actualizar la señal en volts
        function obj = setVolts(obj, volts)
            obj.volts = volts;
        end
      
        % Obtiene la señal en volts
        function getVolts = getVolts(obj)
            getVolts = obj.volts;
        end
        
        %% Señal en gravedades
        % Permite actualizar la señal en gravedades
        function obj = setG(obj, g)
            obj.g = g;
        end
      
        % Obtiene la señal en gravedades
        function getG = getG(obj)
            getG = obj.g;
        end
        
        %% Señal en Gales
        % Permite actualizar la señal en gales
        function obj = setGal(obj, gal)
            obj.gal = gal;
        end
        
        % Obtiene la frecuencia de muestreo
        function getGal = getGal(obj)
            getGal = obj.gal;
        end
        
        %% Frecuencia de Muestreo
        % Obtiene la frecuencia de muestreo
        function getSRate = getSRate(obj)
            getSRate = obj.sampleRate;
        end
        
        %permite acutalizar la frecuencia de muestreo
        function obj = setSRate(obj, sampleRate)
            obj.sampleRate = sampleRate;
        end
        
        % establece los maximos y los minimos de la señal
        function [obj]=maxmin(obj,n)
        % maximos de las componentes
        
        switch n
            case 1
                [max_X,I_X]=max(obj.volts(:,1));
                [max_Y,I_Y]=max(obj.volts(:,2));
                [max_Z,I_Z]=max(obj.volts(:,3));

                [min_X,Imin_X]=min(obj.volts(:,1));
                [min_Y,Imin_Y]=min(obj.volts(:,2));
                [min_Z,Imin_Z]=min(obj.volts(:,3));
            case 2
                [max_X,I_X]=max(obj.g(:,1));
                [max_Y,I_Y]=max(obj.g(:,2));
                [max_Z,I_Z]=max(obj.g(:,3));

                [min_X,Imin_X]=min(obj.g(:,1));
                [min_Y,Imin_Y]=min(obj.g(:,2));
                [min_Z,Imin_Z]=min(obj.g(:,3));
            case 3
                [max_X,I_X]=max(obj.gal(:,1));
                [max_Y,I_Y]=max(obj.gal(:,2));
                [max_Z,I_Z]=max(obj.gal(:,3));

                [min_X,Imin_X]=min(obj.gal(:,1));
                [min_Y,Imin_Y]=min(obj.gal(:,2));
                [min_Z,Imin_Z]=min(obj.gal(:,3));
        end
            %componente X
            if abs(min_X)>max_X
                obj.maxs(1,1)=min_X;
                obj.maxs(1,2)=Imin_X;
            else
                obj.maxs(1,1)=max_X;
                obj.maxs(1,2)=I_X;
            end

            %componente Y
            if abs(min_Y)>max_Y
                obj.maxs(2,1)=min_Y;
                obj.maxs(2,2)=Imin_Y;
            else
                obj.maxs(2,1)=max_Y;
                obj.maxs(2,2)=I_Y;
            end

            %componente Z
            if abs(min_Z)>max_Z
                obj.maxs(3,1)=min_Z;
                obj.maxs(3,2)=Imin_Z;
            else
                obj.maxs(3,1)=max_Z;
                obj.maxs(3,2)=I_Z;
            end    
        end
        
        % Filtra la señal mediante un filtro paso-bajo y otro paso alto de
        % butterworth con la frecuencias y orden definidos
        function [obj] = filtrar(obj, fmin,fmax,n)
            %filtros componente X
            % filtro paso bajo
            Wn=fmax/(obj.sampleRate/2);               %dudaa!!! 
            n=4;
            [b,a] = butter(n,Wn,'low');
            filt_auxX=filtfilt(b,a,obj.gal(:,1));
            
            % filtro paso alto
            Wn=fmin/(obj.sampleRate/2);
            n=4;
            [b,a] = butter(n,Wn,'high');
            filt_X=filtfilt(b,a,filt_auxX);
            
            %filtros componente y
            % filtro paso bajo
            Wn=fmax/(obj.sampleRate/2);               %dudaa!!! 
            n=4;
            [b,a] = butter(n,Wn,'low');
            filt_auxY=filtfilt(b,a,obj.gal(:,2));
            
            % filtro paso alto
            Wn=fmin/(obj.sampleRate/2);
            n=4;
            [b,a] = butter(n,Wn,'high');
            filt_Y=filtfilt(b,a,filt_auxY);
            
            %filtros componente z
            % filtro paso bajo
            Wn=fmax/(obj.sampleRate/2);               %dudaa!!! 
            n=4;
            [b,a] = butter(n,Wn,'low');
            filt_auxZ=filtfilt(b,a,obj.gal(:,3));
            
            % filtro paso alto
            Wn=fmin/(obj.sampleRate/2);
            n=4;
            [b,a] = butter(n,Wn,'high');
            filt_Z=filtfilt(b,a,filt_auxZ);
            
            obj.gal(:,1)=filt_X;
            obj.gal(:,2)=filt_Y;
            obj.gal(:,3)=filt_Z;            
        end
        
        % Método para calcular la frecuencia de fourier
        function [obj] = espectroFourier(obj)
            X=obj.gal(:,1);
            Y=obj.gal(:,1);
            Z=obj.gal(:,1);            
            
            L1=length(X);
            NFFT1 = 2^nextpow2(L1); 
            X_fourier = fft(X,NFFT1)/L1;
            frec = obj.sampleRate/2*linspace(0,1,NFFT1/2+1);
            fourier=figure('Name','Espectro de Fourier','Position',[10 10 650 250]);
            subplot(1,3,1)
%             subplot(3,1,1)
            semilogx(frec,100*abs(X_fourier(1:NFFT1/2+1))) 
            % plot(f1,100*abs(Y1(1:NFFT1/8+1))) 
            title('N-S')
            ylabel('Amplitud')
            axis square;
                        
%             L1=length(Y);
%             NFFT1 = 2^nextpow2(L1); 
            Y_fourier= fft(Y,NFFT1)/L1;
            subplot(1,3,2)
%             subplot(3,1,2)
            semilogx(frec,100*abs(Y_fourier(1:NFFT1/2+1))) 
            % plot(f1,100*abs(Y1(1:NFFT1/8+1))) 
            title('N-S')
            ylabel('Amplitud')
            axis square;
            
%             L1=length(Z);
%             NFFT1 = 2^nextpow2(L1); 
            Z_fourier = fft(Z,NFFT1)/L1;
            subplot(1,3,3)
%             subplot(3,1,3)
            semilogx(frec,100*abs(Z_fourier(1:NFFT1/2+1))) 
            % plot(f1,100*abs(Y1(1:NFFT1/8+1))) 
            title('N-S')
            ylabel('Amplitud')
            axis square;
            
            obj.fourier(:,1)=X_fourier;
            obj.fourier(:,2)=Y_fourier;
            obj.fourier(:,3)=Z_fourier;
            obj.fourier_HZ=frec;
            
%             F = getframe(fourier);
%             imwrite(F.cdata,'Fourier.jpg');
        end
        
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
    end

end