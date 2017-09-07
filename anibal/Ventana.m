function []=Ventana()
prompt = {'Capital:','Estación:','Ipocentro:','Magnitud:','Profundidad:',...
    'Fecha:','Hora:','Geología:','Topografía:'};
title = 'Generar Reporte';
lines = 1;
% def = {'my_image','hsv'};
answer = inputdlg(prompt,title,lines);
% capital=answer{1};
% estacion=answer{2};
% ipocentro=answer{3};
% magnitud=answer{4};
% profundidad=answer{5};
% fecha=answer{6};
% hora=answer{7};
% geologia=answer{8};
% topo=answer{9};

assignin('base','capital',answer{1});
assignin('base','estacion',answer{2});
assignin('base','ipocentro',answer{3});
assignin('base','magnitud',answer{4});
assignin('base','profundidad',answer{5});
assignin('base','fecha',answer{6});
assignin('base','hora',answer{7});
assignin('base','geologia',answer{8});
assignin('base','topo',answer{9});

