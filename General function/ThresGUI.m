function varargout = ThresGUI(varargin)
%
% THRESGUI
%
% THRESGUI IS A GUI CREATED TO VISUALIZATION OF THRESHOLDING PROCESS OF AN
% IMAGE. YOU CAN CHOICE THE THRESHOLDING METHOD OR ARBITRARY THRESHOLDS.
% THE THRESHOLDING METHODS AVALIABLE ARE:
% .GRAYTRHESH (MATLAB)
% .KAPUR
% .ITERATIVE
% .TRIANGULAR
%
% REQUIREMENTS: MATLAB 7.0, RELEASE 14.
%
% CREATED BY:
% IVAN DARIO ORDOÑEZ SEPULVEDA
% CHEMICAL ENGINEER, 2006
%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThresGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ThresGUI_OutputFcn, ...
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

% --- Executes just before ThresGUI is made visible.
function ThresGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%========================================================================
%========================================================================

% --- Outputs from this function are returned to the command line.
function varargout = ThresGUI_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
global Archivo NameArchivo
%[NameArchivo,PathArchivo] = uigetfile({'*.jpg';, '*.tif';, '*.gif';, '*.bmp';, '*.png';, ...
%    '*.hdf';, '*.pcx';, '*.xwd';, '*.ico';, '*.cur';, '*.ras';, ...
%    '*.pbm';, '*.pgm';, '*.ppm'});
set(handles.Save,'enable','off')
set(handles.TipoSegmentacion,'enable','off')
[NameArchivo,PathArchivo] = uigetfile({'*.jpg;*.tif;*.gif;*.bmp;*.png; *.hdf;*.pcx;*.xwd;*.ico;*.cur;*.ras;*.pbm;*.pgm;*.ppm'});
if ~isequal(NameArchivo, 0)
    Archivo=[PathArchivo,NameArchivo];
    CargarImagen(handles)
    set(handles.Save,'enable','on')
    set(handles.TipoSegmentacion,'enable','on')
    set(handles.DefinidoUsuario,'checked','on')
end

% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
global ImagenUmbral
if isempty(ImagenUmbral)==1,msgbox('Doesn''t exist an image');return,end
[filename,pathname] = uiputfile({'*.jpg';, '*.tif';, '*.gif';, '*.bmp';, '*.png';, ...
    '*.hdf';, '*.pcx';, '*.xwd';, '*.ico';, '*.cur';, '*.ras';, ...
    '*.pbm';, '*.pgm';, '*.ppm'},'Save file name');
if isequal(filename,0) | isequal(pathname,0)
   errordlg('Saving canceled','Threshold GUI'); error('Saving canceled')
else
    try
        imwrite(ImagenUmbral,[pathname,filename]);
    catch
        errordlg('Error during saving','Threshold GUI'),error('Error during saving')
    end %try
end %if

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.figure1)

% --------------------------------------------------------------------
function graythresh_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
global ImagenGris 
NoChecked(handles)
set(handles.graythresh,'checked','on')
Umbral = graythresh(ImagenGris)*255;
GraficarHistoUmbral(Umbral,handles,'all')

% --------------------------------------------------------------------
function kapur_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
global ImagenGris 
NoChecked(handles)
set(handles.kapur,'checked','on')
Umbral = Kapur1(ImagenGris)*255;
GraficarHistoUmbral(Umbral,handles,'all')

% --------------------------------------------------------------------
function triangular_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
global ImagenGris 
NoChecked(handles)
set(handles.triangular,'checked','on')
Umbrales = Triangular1(ImagenGris)*255;
set(handles.triangular,'userdata',Umbrales)
GraficarHistoUmbral(Umbrales(1),handles,'trian')

% --------------------------------------------------------------------
function iterativo_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
global ImagenGris 
NoChecked(handles)
set(handles.iterativo,'checked','on')
Umbral = Iterativo1(ImagenGris)*255;
GraficarHistoUmbral(Umbral,handles,'all')

% --------------------------------------------------------------------
function DefinidoUsuario_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
NoChecked(handles)
set(handles.DefinidoUsuario,'checked','on')
set(handles.Umbralizar,'enable','on')
set(handles.Umbralizar,'string','Thresholding')
if get(handles.ActualAutomatico,'value')==1
    set(handles.Umbralizar,'enable','off')
end %if

% --- Executes on button press in Umbralizar.
%-------------------------------------------------------------------------
function Umbralizar_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
if isequal(get(handles.triangular,'checked'),'off')==1
    Umbral=get(handles.slider1,'value');
    GraficarUmbral(Umbral,handles);
    NoChecked(handles);set(handles.DefinidoUsuario,'checked','on')
else
    UmbralActual=get(handles.slider1,'value');
    Umbrales=get(handles.triangular,'userdata');
    Pos=find(Umbrales==UmbralActual);
    if isempty(Pos),msgbox('Error de Ejecución');return,end
    if Pos==length(Umbrales),Pos=0;,end %if
    GraficarHistoUmbral(Umbrales(Pos+1),handles,'')
end %if

% --- Executes on slider movement.
%-------------------------------------------------------------------------
function slider1_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
set(handles.Intensidad,'visible','on')
Umbral=get(hObject,'value');
GraficarHistograma(Umbral,handles)
if get(handles.ActualAutomatico,'value')==1
    GraficarUmbral(Umbral,handles);
end %if
NoChecked(handles);set(handles.DefinidoUsuario,'checked','on')

% --- Executes during object creation, after setting all properties.
%-------------------------------------------------------------------------
function slider1_CreateFcn(hObject, eventdata, handles)
%-------------------------------------------------------------------------
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
set(hObject,'SliderStep',[1/255 1/255]);%Saltra cada 1

% --- Executes on button press in ActualAutomatico.
%-------------------------------------------------------------------------
function ActualAutomatico_Callback(hObject, eventdata, handles)
%-------------------------------------------------------------------------
set(handles.Umbralizar,'enable','on')
set(handles.Umbralizar,'string','Thresholding')
set(handles.Umbralizar,'enable','off')
if get(hObject,'Value')==1
    set(handles.Umbralizar,'enable','off')
else
    set(handles.Umbralizar,'enable','on')
end %if
NoChecked(handles);set(handles.DefinidoUsuario,'checked','on')

% --------------------------------------------------------------------
function TipoSegmentacion_Callback(hObject, eventdata, handles)


%==========================================================================
function CargarImagen(handles)
%==========================================================================
global Archivo NameArchivo ImagenGris

subplot(221);img=imread(Archivo);imshow(img);title('Original Image');axis off
try
    if Esgray(img)==0
        imggris=rgb2gray(img);subplot(222);imshow(imggris);title('Gray Scale Image')
    else
        imggris=img;subplot(222);imshow(imggris);title('Gray Scale Image')
    end %if
catch
    errordlg('Error during image processing','Threshold GUI');error('Error during image processing')
end %try
subplot(223);imhist(imggris);title('Histogram')
ImagenGris=imggris;

Min=0;Max=255;
set(handles.slider1,'visible','on')
set(handles.slider1,'Min',Min)
set(handles.slider1,'Max',Max)
set(handles.slider1,'value',Min)
set(handles.Umbralizar,'visible','on')
set(handles.ActualAutomatico,'visible','on')
set(handles.text2,'visible','on')
set(handles.Intensidad,'string',int2str(Min))
set(handles.Archivo,'visible','on')
set(handles.Archivo,'string',NameArchivo)
GraficarUmbral(Min,handles)


%==========================================================================
function GraficarUmbral(Umbral,handles)
%==========================================================================
global ImagenGris ImagenUmbral
imgUmbral=im2bw(ImagenGris,Umbral/255);
subplot(224);
imshow(imgUmbral)
title(['Image segmented in threshold = ',int2str(Umbral)])
ImagenUmbral=imgUmbral;

%==========================================================================
function GraficarHistograma(Umbral,handles)
%==========================================================================
global ImagenGris
persistent H
set(handles.Intensidad,'String',Umbral)
subplot(223);imhist(ImagenGris);title('Histogram')
Escala=axis;
if ishandle(H)==1,delete(H),end%if
H=line([Umbral Umbral],[Escala(3:4)],'color','r');

%==========================================================================
function GraficarHistoUmbral(Umbral,handles,Tipo)
%==========================================================================
GraficarHistograma(Umbral,handles)
GraficarUmbral(Umbral,handles)
set(handles.slider1,'value',Umbral)
switch Tipo
    case 'all' %GRAYTHRESH ITERATIVO KAPUR
        set(handles.ActualAutomatico,'value',0)
        set(handles.Umbralizar,'enable','on')
        set(handles.Umbralizar,'string','Thresholding')
        set(handles.Umbralizar,'enable','off')
    case 'trian' %TRIANGULAR
        set(handles.ActualAutomatico,'value',0)
        set(handles.Umbralizar,'enable','on')
        set(handles.Umbralizar,'string','next Threshold >>')
end %switch

%==========================================================================
function NoChecked(handles)
%==========================================================================
set(handles.graythresh,'checked','off')
set(handles.kapur,'checked','off')
set(handles.triangular,'checked','off')
set(handles.iterativo,'checked','off')
set(handles.DefinidoUsuario,'checked','off')


%==========================================================================
function umbral = Kapur1(Imagen)
%==========================================================================
%AUTOR OF THIS SUBROUTINE: BIOMEDICAL ENGINEERING RESEARCH GROUP - UIS
if ~Esgray(Imagen)                                  %Validar si una Imagen se encuentra en Escala de Grises o en color
    uiwait(msgbox('La Imagen no se encuentra en escala de grises','Error','modal'))
    umbral = 0;
else 
    [fil col] = size(Imagen);
    Pixeles = fil * col;
    h1 = imhist(Imagen);
    Pi = h1/Pixeles;                            %Calcular las probabilidades de cada escala de gris
    Pt = zeros(256,1);
    Pt(1) = Pi(1);
    for i = 2:256                               %Calcular la Función de Probabilidad Acumulada
        Pt(i)=Pt(i-1)+Pi(i);
    end

    Hb = zeros(1,256);                          %Calcular Hb(u) y Hw(u) donde la suma de los dos es máxima suponiendo un umbral
    Hw = zeros(1,256);
    for i = 1:256
        if Pt(i) > 0
            for j = 1 : i
                if Pi(j) > 0
                    Hb(i) = Hb(i) + ((Pi(j) / Pt(i)) * log(Pi(j) / Pt(i)));
                end
            end
        end
    end

    for i = 1:256
        if (1-Pt(i)) > 0
            for j = i + 1 : 256
                if Pi(j) > 0
                    Hw(i) = Hw(i) + ((Pi(j) / (1-Pt(i))) * log(Pi(j) / (1-Pt(i))));
                end
            end
        end
    end
    
    Hb = -Hb;
    Hw = -Hw;
    H = Hb + Hw;
    [a, b] = max(H(:));
    umbral = b-1;
    umbral = umbral/255;
end

%==========================================================================
function Umbral = Triangular1(Imagen);
%==========================================================================
%AUTOR OF THIS SUBROUTINE: BIOMEDICAL ENGINEERING RESEARCH GROUP - UIS
if ~Esgray(Imagen)                                      %Validar si una Imagen se encuentra en Escala de Grises o en color
    uitwait(msgbox('La imagen no se encuentra en escala de grises','Error','modal'))
    Umbral = 0;
else
    z = 0;
    H = imhist(Imagen);
    H = Promedio(H);
    dh = Derivada(H);
    [maxim minim]= MaxiMini(dh);
    for i = 1 : length(maxim)-1;
        clear D, k = 0;
        x1 = maxim(i); y1 = H(maxim(i));
        x2 = maxim(i + 1); y2 = H(maxim(i + 1));
        M = (y2 - y1) / (x2 - x1);
        if (x2 - 1) - (x1 + 1) > 0
            for j = x1 : x2;
                Px = j;
                Py = H(j);
                k = k + 1;
                D(k,1) = sqrt(((Px - x1) * M - Py + y1)^2 / (M^2 + 1));
                D(k,2) = j;
            end
            [Pu p]= max(D(:,1));
            if Pu > 0
                z = z + 1;
                Umbral(z) = D(p,2);
            end
        end
    end
    Umbral = Repetidos(Umbral,10)/255;
end

function x = Promedio(h);
n = 5;
g = zeros(256 + 2 * round(n / 2 - 1), 1);
g(1 + round(n / 2 - 1) : length(g) - round(n / 2 - 1)) = h;
for i = 1 + round(n / 2 - 1) : length(g)-round(n / 2 - 1);
    sum = 0;
    for j = i - round(n / 2 - 1) : i + round(n / 2 - 1)
        sum = sum + g(j);
    end
    x(i - round(n / 2 - 1), 1) = round(sum / n);
end

function dx = Derivada(h);
dx=zeros(256,1);
for i=3:254;
    dx(i)= (h(i-2) - 8*(h(i-1)) + 8*(h(i+1)) - h(i+2))/12;
end

function [Maxi, Mini] = MaxiMini(d);
j = 1; y = 1;
for i = 2 : 256;
    if (d(i - 1) < 0 && d(i) >= 0)
        Mini(j) = i;
        j = j + 1;
    elseif (d(i - 1) >= 0 && d(i) < 0)
        Maxi(y) = i;
        y = y + 1;
    end
end

function Vector = Repetidos(Arreglo, D);
N_Ind = length(Arreglo);
for i = 1 : N_Ind - 1;
    if Arreglo(i) > 0
        for j = i + 1 : N_Ind;
            if Arreglo(i) == Arreglo(j) || (abs(Arreglo(i)-Arreglo(j))<D)
                Arreglo_Nuevo(j) = 0;
            else
                Arreglo_Nuevo(j) = Arreglo(j);
            end
        end
    end
end
[i j Vector] = find(Arreglo_Nuevo);

%==========================================================================
function Umbral = Iterativo1(Imagen)
%==========================================================================
%AUTOR OF THIS SUBROUTINE: BIOMEDICAL ENGINEERING RESEARCH GROUP - UIS
if ~Esgray(Imagen)                                      %Validar si una Imagen se encuentra en Escala de Grises o en color
    uitwait(msgbox('La imagen no se encuentra en escala de grises','Error','modal'))
    Umbral = 0;
else 
    Histograma = imhist(Imagen);
    Grises = find(Histograma);                      %Localizar los niveles diferentes a cero
    Maximo = max(Grises);                           %Nivel máximo de gris
    Minimo = min(Grises);                           %Nivel mínimo de gris
    Umbral = Minimo + (Maximo - Minimo) / 2;     %Umbral Inicial
    Umbralp = 0;                                    %Umbral Anterior
    while (abs(Umbral - Umbralp) > 1)
        Mayores = find(Imagen > Umbral);
        Menores = find(Imagen <= Umbral);
        m1 = mean(Imagen(Mayores));
        m2 = mean(Imagen(Menores));
        Umbralp = Umbral;
        Umbral = round((m1+m2)/2);              %Calcular el umbral óptimo
    end
    Umbral = Umbral/255;
end

%==========================================================================
function y = Esgray(x)
%==========================================================================
y = ndims(x)==2 && ~isempty(x);

if islogical(x)
  y = false;
elseif ~isa(x, 'uint8') && ~isa(x, 'uint16') && y
  % At first just test a small chunk to get a possible quick negative
  [m,n] = size(x);
  chunk = x(1:min(m,10),1:min(n,10));         
  y = min(chunk(:))>=0 && max(chunk(:))<=1;
  % If the chunk is an intensity image, test the whole image
  if y
    y = min(x(:))>=0 && max(x(:))<=1;
  end
end 

