function varargout = RAC_GUI(varargin)
% RAC_GUI M-file for RAC_GUI.fig
%      RAC_GUI, by itself, creates a new RAC_GUI or raises the existing
%      singleton*.
%
%      H = RAC_GUI returns the handle to a new RAC_GUI or the handle to
%      the existing singleton*.
%
%      RAC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAC_GUI.M with the given input arguments.
%
%      RAC_GUI('Property','Value',...) creates a new RAC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RAC_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RAC_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RAC_GUI

% Last Modified by GUIDE v2.5 10-May-2013 17:40:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RAC_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RAC_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before RAC_GUI is made visible.
function RAC_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RAC_GUI (see VARARGIN)

% Choose default command line output for RAC_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RAC_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RAC_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_Cargar.
function button_Cargar_Callback(hObject, eventdata, handles)
% hObject    handle to button_Cargar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global station
%% selecciona archivo a cargar
[nom_arch, ruta]=uigetfile({'*.EVT';'*.*'}, 'Seleccione el evento a procesar:');
ejecutar='!kw2asc32.exe'
comand=strcat(ejecutar,{' '},nom_arch);
% comand=strcat(ejecutar,{' '},ruta,nom_arch);
eval(comand{1});

%% guarda y corrigue la linea base de las trazas 
% componente X
acelero(:,1)=load(strcat(ruta,nom_arch(1:end-4),'.001'));
% componente Y
acelero(:,2)=load(strcat(ruta,nom_arch(1:end-4),'.002'));
% componente Z
acelero(:,3)=load(strcat(ruta,nom_arch(1:end-4),'.003'));

%% variables dependientes de la configuracion del equipo
%sensibilidad
gen_k=1.25;
%constante para conversion de gravedades a gales
gal_k=981;
%frecuencia de muestreo
fm=200;

%% crea el objeto con las trazas 
station = Acelerograma(acelero, gen_k, gal_k,fm);

%% crea la imagen inicial para guardarla

%% graficar objeto dependiendo de la vista seleccionada
val = get(handles.popupmenu_view,'Value');
station= maxmin(station,val);
switch val
    case 1
%         axes(handles.axesX);        
        subplot(3,1,1)
        title('Acelerograma');  
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.volts(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2)
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.volts(:,2),'MarkerSize',17)
        ylabel ('Aceleración (Volts)')
        title ('E-O') 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.volts(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
        
    case 2
        
%         axes(handles.axesX);
        subplot(3,1,1)
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.g(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2)
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.g(:,2),'MarkerSize',17)
        ylabel ('Aceleración (m/seg/seg)')
        title ('E-O') 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.g(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
        
    case 3
%         axes(handles.axesX);       
        subplot(3,1,1)
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.gal(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2)
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.gal(:,2),'MarkerSize',17)
        ylabel ('Aceleración (cm/seg/seg)')
        title ('E-O') 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.gal(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
        
end

% --- Executes on selection change in popupmenu_view.
function popupmenu_view_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_view contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_view
global station;
val = get(hObject,'Value');
station= maxmin(station,val);
 switch val 
     case 1
%         axes(handles.axesX);
        subplot(3,1,1);         
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.volts(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2);
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.volts(:,2),'MarkerSize',17)
        ylabel ('Aceleración (Volts)')
        title ('E-O') 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.volts(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
     
    case 2
%         axes(handles.axesX);
        subplot(3,1,1)
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.g(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2)
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.g(:,2),'MarkerSize',17)
        ylabel ('Aceleración (m/seg/seg)')
        title ('E-O') 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.g(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
        
    case 3
%         axes(handles.axesX); 
        subplot(3,1,1);
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.gal(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2);
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.gal(:,2),'MarkerSize',17)
        ylabel ('Aceleración (cm/seg/seg)')
        title ('E-O') 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3);
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.gal(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
        
end

% --- Executes during object creation, after setting all properties.
function popupmenu_view_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_process.
function popupmenu_process_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_process contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_process
global station;
val = get(hObject,'Value');

switch val
    %% Filtrado
    case 2       
        % variables para el filtrado
        f_minimo = str2double(get(handles.edit_fmin,'string'));        
        f_maximo=str2double(get(handles.edit_fmax,'string'));
        orden=str2double(get(handles.edit_orden,'string'));

        % filtrar señales
        station = filtrar(station, f_minimo,f_maximo,orden);
        station= maxmin(station,3);        
                
        %% graficar en la GUI
%         axes(handles.axesX); 
        subplot(3,1,1)
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.gal(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2)
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.gal(:,2),'MarkerSize',17)
        ylabel ('Aceleración (cm/seg/seg)')
        title ('E-O') 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.gal(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
    
    %% calcular espectro de fourier
    case 3        
        [station] = espectroFourier(station);
        
    case 4
        %% espectro de respuesta
        % variables
        ro=0.05; % Amortiguacion
        T_pas=0:0.01:4; % Paso para el periodo
        wn_pas=2*pi./T_pas; % Frecuencia

        [station] = espectroRespuesta(station,ro,wn_pas,T_pas);    

end

% --- Executes during object creation, after setting all properties.
function popupmenu_process_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_reporte.
function button_reporte_Callback(hObject, eventdata, handles)
% hObject    handle to button_reporte (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global station
val = get(handles.popupmenu_process,'Value');

if(val<2)
    % variables para el filtrado
    f_minimo = str2double(get(handles.edit_fmin,'string'));        
    f_maximo=str2double(get(handles.edit_fmax,'string'));
    orden=str2double(get(handles.edit_orden,'string'));

    % filtrar señales
    station = filtrar(station, f_minimo,f_maximo,orden);
    station= maxmin(station,3);
end

if(val<3)
    %% Fourier
    [station] = espectroFourier(station);
    
end

if(val<4)
    %% espectro de respuesta
    % variables
    ro=0.05; % Amortiguacion
    T_pas=0:0.01:4; % Paso para el periodo
    wn_pas=2*pi./T_pas; % Frecuencia    
    [station] = espectroRespuesta(station,ro,wn_pas,T_pas);    
end   

% close all;
Ventana()
assignin('base','station_report',station);
assignin('base','Tp',T_pas);
report RAC


function edit_fmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fmax as text
%        str2double(get(hObject,'String')) returns contents of edit_fmax as a double


% --- Executes during object creation, after setting all properties.
function edit_fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fmin as text
%        str2double(get(hObject,'String')) returns contents of edit_fmin as a double


% --- Executes during object creation, after setting all properties.
function edit_fmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_orden_Callback(hObject, eventdata, handles)
% hObject    handle to edit_orden (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_orden as text
%        str2double(get(hObject,'String')) returns contents of edit_orden as a double


% --- Executes during object creation, after setting all properties.
function edit_orden_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_orden (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Menu_Archivo_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Archivo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Cargar_EVT_Callback(hObject, eventdata, handles)
% hObject    handle to Cargar_EVT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global station
%% selecciona archivo a cargar
[nom_arch, ruta]=uigetfile({'*.EVT';'*.*'}, 'Seleccione el evento a procesar:');
ejecutar='!kw2asc32.exe'
comand=strcat(ejecutar,{' '},nom_arch);
% comand=strcat(ejecutar,{' '},ruta,nom_arch);
eval(comand{1});

%% guarda y corrigue la linea base de las trazas 
% componente X
acelero(:,1)=load(strcat(ruta,nom_arch(1:end-4),'.001'));
% componente Y
acelero(:,2)=load(strcat(ruta,nom_arch(1:end-4),'.002'));
% componente Z
acelero(:,3)=load(strcat(ruta,nom_arch(1:end-4),'.003'));

%% variables dependientes de la configuracion del equipo
%sensibilidad
gen_k=1.25;
%constante para conversion de gravedades a gales
gal_k=981;
%frecuencia de muestreo
fm=200;

%% crea el objeto con las trazas 
station = Acelerograma(acelero, gen_k, gal_k,fm);

%% crea la imagen inicial para guardarla

%% graficar objeto dependiendo de la vista seleccionada
val = get(handles.popupmenu_view,'Value');
station= maxmin(station,val);
switch val
    case 1
%         axes(handles.axesX);        
        subplot(3,1,1)
        title('Acelerograma');  
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.volts(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2)
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.volts(:,2),'MarkerSize',17)
        ylabel ('Aceleración (Volts)')
        title ('E-O') 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.volts(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (Volts):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
        
    case 2
        
%         axes(handles.axesX);
        subplot(3,1,1)
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.g(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2)
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.g(:,2),'MarkerSize',17)
        ylabel ('Aceleración (m/seg/seg)')
        title ('E-O') 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.g(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (m/seg/seg):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
        
    case 3
%         axes(handles.axesX);       
        subplot(3,1,1)
        plot(station.timeTotal(station.maxs(1,2)),station.maxs(1,1),'.r',station.timeTotal,station.gal(:,1),'MarkerSize',17)
        title ('N-S'); 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(1,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,2)
        plot(station.timeTotal(station.maxs(2,2)),station.maxs(2,1),'.r',station.timeTotal,station.gal(:,2),'MarkerSize',17)
        ylabel ('Aceleración (cm/seg/seg)')
        title ('E-O') 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(2,1))));
        set(h,'Interpreter','none')
        
        subplot(3,1,3)
        plot(station.timeTotal(station.maxs(3,2)),station.maxs(3,1),'.r',station.timeTotal,station.gal(:,3),'MarkerSize',17)
        xlabel ('Tiempo (seg)')  
        title ('Ver') 
        h = legend(strcat('MAX (cm/seg/seg):',{' '},num2str(station.maxs(3,1))));
        set(h,'Interpreter','none') 
        
end

% --------------------------------------------------------------------
function Cargar_Seisan_Callback(hObject, eventdata, handles)
% hObject    handle to Cargar_Seisan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% selecciona archivo a cargar
% componente X
[nom_arch, ruta]=uigetfile({'*.EVT';'*.*'}, 'Seleccione la componente Vertical (Z) del evento:');
acelero_X=dlmread(strcat(ruta,nom_arch));

% componente Y
[nom_arch, ruta]=uigetfile({'*.EVT';'*.*'}, 'Seleccione la componente Norte-Sur (X) del evento:');
acelero_Y=dlmread(strcat(ruta,nom_arch));

% componente Z
[nom_arch, ruta]=uigetfile({'*.EVT';'*.*'}, 'Seleccione la componente Este-Oeste (Y) del evento:');
acelero_Z=dlmread(strcat(ruta,nom_arch));
