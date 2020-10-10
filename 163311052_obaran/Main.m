
function varargout = Main(varargin)
% MAIN MATLAB code for Main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main

% Last Modified by GUIDE v2.5 23-Jun-2020 20:00:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_OutputFcn, ...
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


% --- Executes just before Main is made visible.
function Main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main (see VARARGIN)

% Choose default command line output for Main
handles.output = hObject;
axes(handles.axes3)
imshow('mylogo.png');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function dislem_Callback(hObject, eventdata, handles)
% hObject    handle to dislem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function nbaz_Callback(hObject, eventdata, handles)
% hObject    handle to nbaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fislem_Callback(hObject, eventdata, handles)
% hObject    handle to fislem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mislem_Callback(hObject, eventdata, handles)
% hObject    handle to mislem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gsislem_Callback(hObject, eventdata, handles)
% hObject    handle to gsislem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function asinma_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
A=image;
A=im2bw(A);
B=[1 1 0];
C=padarray(A,[0 1],1);
D=false(size(A));
for i=1:size(C,1)
    for j=1:size(C,2)-2
        L=C(i,j:j+2);
        K=find(B==1);
       if(L(K)==1)
        D(i,j)=1;
        end
    end
end
axes(handles.axes2)
imshow(D);


% --------------------------------------------------------------------
function genisleme_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
A=image;
A=im2bw(A);
B=[1 1 1 1 1 1 1;];
C=padarray(A,[0 3]);
D=false(size(A));
for i=1:size(C,1)
    for j=1:size(C,2)-6
        D(i,j)=sum(B&C(i,j:j+6));
    end
end
axes(handles.axes2)
imshow(D);

% --------------------------------------------------------------------
function acma_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
I = image;
imshow(I);
se = strel('disk',5);
I_opened = imopen(I,se);
axes(handles.axes2)
imshow (I_opened);

% --------------------------------------------------------------------
function kapama_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
originalBW = image;
se = strel('disk',10);
closeBW = imclose(originalBW,se);
axes(handles.axes2)
imshow(closeBW)


% --------------------------------------------------------------------
function alcakgeciren_Callback(hObject, eventdata, handles)
% hObject    handle to alcakgeciren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function yuksekgeciren_Callback(hObject, eventdata, handles)
% hObject    handle to yuksekgeciren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function esik_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',true);
set(handles.slider1, 'max', 255); %E?ik de?eri minimum 0 maksimum 255 olabilir.
global image
x=image;
x=rgb2gray(x);
x=double(x);
tot=0;
[a,b]=size(x);
y=zeros(a,b);
for i=1:a
    for j=1:b
        y(i,j)=0;
    end
end
for i=1:a
    for j=1:b
        tot=tot+x(i,j);
    end
end
sliderValue = get(handles.slider1, 'Value');
thr=sliderValue; % Slider de?erini e?ik de?eri olarak ayarl?yoruz.
for i=1:a
    for j=1:b
        if x(i,j) > thr %Bu k?s?mda piksel de?eri verdi?imiz e?ik de?erinden yüksekse 
            y(i,j)=0;   %piksel de?erini 0 dü?ük veya e?itse piksel de?erini 1 yap?yor.
        else
            y(i,j)=1;
        end
    end
end
axes(handles.axes2)
imshow(y);



% --------------------------------------------------------------------
function ldon_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',true);
set(handles.slider1,'max',5);%Burada slider maximum de?eri 5 olmu?tur.
deger=get(handles.slider1,'Value');
if deger>5
    set(handles.slider1,'Value',3);
end
global image
a1 = image;
a = double(a1)/255; 
c = deger;
f = c*log(1 + (a)); % Log Dönü?ümü
axes(handles.axes2);       
imshow(f);

% --------------------------------------------------------------------
function gsdon_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',true);
deger=get(handles.slider1,'Value');
if deger>5
    set(handles.slider1,'Value',5);
    set(handles.slider1,'max',5);
    deger=1;
end
global image;
F=rgb2gray(image);
C=0.01;
gamma=deger;
I=uint8(C.*((double(F)).^gamma));
axes(handles.axes2);
imshow(I);
  

% --------------------------------------------------------------------
function kgerme_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
itemp = image; 
i = itemp(:,:,1);
rtemp = min(i);     
rmin = min(rtemp);
rtemp = max(i); 
rmax = max(rtemp);
m = 255/(rmax - rmin);
c = 255 - m*rmax; 
i_new = m*i + c; 
axes(handles.axes2)
imshow(i_new);

% --------------------------------------------------------------------
function hesitleme_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
orginal =image;
[rows,columns,~] = size(orginal);
finalResult = uint8(zeros(rows,columns));
pixelNumber = rows*columns;
frequncy = zeros(256,1);
pdf = zeros(256,1);
cdf = zeros(256,1);
cummlative = zeros(256,1);
outpic = zeros(256,1);
for i = 1:1:rows
    for j = 1:1:columns
        val = orginal(i,j);
        frequncy(val+1) = frequncy(val+1)+1;
        pdf(val+1) = frequncy(val+1)/pixelNumber;
    end
end
sum =0 ;
intensityLevel = 255;
for i = 1:1:size(pdf)
    sum =sum +frequncy(i);
    cummlative(i) = sum;
    cdf(i) = cummlative(i)/ pixelNumber;
    outpic(i) = round(cdf(i) * intensityLevel);
end
for i = 1:1:rows
    for j = 1:1:columns
        finalResult(i,j) = outpic(orginal(i,j) + 1);
    end
end
axes(handles.axes2)
imshow(finalResult);



% --------------------------------------------------------------------
function yresim_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
global filename
global pathname
[filename, pathname]=uigetfile({'*.jpg';'*.gif';'*.png';'*.ppm';'*.pgm'}, 'Görüntü Dosyasini Seç:');
handles.filename=filename;
handles.pathname=pathname;
if filename==0
    msgbox(sprintf('Lütfen bir resim seçiniz.'),'HATA','Error');
else
set(handles.kaydet,'Visible',true);
set(handles.nbaz,'Visible',true);
set(handles.fislem,'Visible',true);
set(handles.mislem,'Visible',true);
set(handles.gsislem,'Visible',true);
axes(handles.axes1)
image=imread([pathname,filename]);
imshow(image);
set(handles.resimadi,'String',filename);
set(handles.resimuzantisi,'String',pathname);
end
% --------------------------------------------------------------------
function kaydet_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
imsave;

% --------------------------------------------------------------------
function sobel_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
A=image;
B=rgb2gray(A);
C=double(B);
for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        Gx=((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        Gy=((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
        B(i,j)=sqrt(Gx.^2+Gy.^2);
    end
end
axes(handles.axes2);
imshow(B);

% --------------------------------------------------------------------
function prewitt_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
openImage = rgb2gray(image);
openImage = im2double(openImage);
b=[-1 -1 -1;0 0 0;1 1 1]/6;
c=[-1 0 1; -1 0 1; -1 0 1]/6;
Gx=abs(conv2(openImage,c,'same'));
Gy=abs(conv2(openImage,b,'same'));
G = sqrt( Gx.^2 + Gy.^2);
out = G > 0.08995; 
axes(handles.axes2);
imshow(out);

% --------------------------------------------------------------------
function laplace_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
a=image;
a=rgb2gray(a);
[r c]=size(a)
a=im2double(a);
filter=[-1 -1 -1;-1 8 -1; -1 -1 -1];
result=a;
for i=2:r-1
    for j=2:c-1
        sum=0;
        row=0;
        col=1;
        
        for k=i-1:i+1
            row=row+1;
            col=1;
            for l=j-1:j+1
                sum = sum+a(k,l)*filter(row,col);               
                col=col+1;
            end
        end
      result(i,j)=sum;      
    end
end
axes(handles.axes2);
imshow(result);
% --------------------------------------------------------------------
function roberts_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
input_image = image; 
input_image = uint8(input_image);  
input_image = rgb2gray(input_image); 
input_image = double(input_image); 
filtered_image = zeros(size(input_image)); 
Mx = [1 0; 0 -1]; 
My = [0 1; -1 0]; 
for i = 1:size(input_image, 1) - 1 
    for j = 1:size(input_image, 2) - 1 
        Gx = sum(sum(Mx.*input_image(i:i+1, j:j+1))); 
        Gy = sum(sum(My.*input_image(i:i+1, j:j+1))); 
        filtered_image(i, j) = sqrt(Gx.^2 + Gy.^2); 
    end
end
filtered_image = uint8(filtered_image); 
axes(handles.axes2);
imshow(filtered_image);

% --------------------------------------------------------------------
function mean_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
I= image;
I = rgb2gray(I);
[R,C] = size(I);
SI = uint8(zeros(R,C));
for r=1:R
for c=1:C
ListOfValues = 1:9;% Çekirdek matris 3*3=9 birimden olu?tu?u için 1/9 kullan?yoruz.
counter = 0;
for i=-1:1
for j=-1:1
counter = counter + 1;
if(((r+i)>0) && ((c+j)>0) && ((r+i)<=R) && ((c+j)<=C))
ListOfValues(counter) = I(r+i,c+j);%Pikselin etraf?ndaki de?erleri listemize aktar?yoruz.
end
end
end
SI(r,c) = median(ListOfValues);
end
end
axes(handles.axes2);
imshow(SI);

% --------------------------------------------------------------------
function median_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image
 A = image;
A = im2double(A);
[m n] = size(A);
Med = [];
for i=2:m-1
    for j=2:n-1 %Resim üzerinde gezerek 9 karelik bir matris üzerinde median i?lemi uygular.
            Med(1) = A(i-1,j-1);
            Med(2) =A(i-1,j) ;
            Med(3) = A(i-1,j+1);
            Med(4) = A(i,j-1);
            Med(5) = A(i,j+1);
            Med(6) = A(i+1, j-1);
            Med(7) = A(i+1,j);
            Med(8) = A(i+1,j+1);
            A(i,j) = median(Med);
    end
end 
axes(handles.axes2);
imshow(A);

% --------------------------------------------------------------------
function gauss_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',true);
deger=get(handles.slider1,'Value');
if deger>50
    set(handles.slider1,'Value',3);
end
set(handles.slider1,'max',50);
global image
I =image;
Iblur = imgaussfilt(I, deger);
axes(handles.axes2);
imshow(Iblur);

function resimadi_Callback(hObject, eventdata, handles)
% hObject    handle to resimadi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resimadi as text
%        str2double(get(hObject,'String')) returns contents of resimadi as a double


% --- Executes during object creation, after setting all properties.
function resimadi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resimadi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resimuzantisi_Callback(hObject, eventdata, handles)
% hObject    handle to resimuzantisi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resimuzantisi as text
%        str2double(get(hObject,'String')) returns contents of resimuzantisi as a double


% --- Executes during object creation, after setting all properties.
function resimuzantisi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resimuzantisi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function svd_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global pathname;
global filename;
global image;
I=image;   
        red = double(I(:,:,1));
       green = double(I(:,:,2));
        blue = double(I(:,:,3)); 
       [u,s,v] = svds(red);
        imred = uint8(u * s * transpose(v)); 
       [u,s,v] = svds(green);
        imgreen = uint8(u * s * transpose(v));
        [u,s,v] = svds(blue);
       imblue = uint8(u * s * transpose(v)); 
     im(:,:,1) = imred;
       im(:,:,2) = imgreen;
       im(:,:,3) = imblue;
    infoEski=imfinfo([pathname,filename]);
imwrite(im,filename,'jpg');
    infoYeni=imfinfo(filename);
    delete(filename);
    set(handles.edit3,'String',"Eski dosya boyutu:"+infoEski.FileSize+" Yeni dosya boyutu:"+infoYeni.FileSize);
 axes(handles.axes2);
 imshow(im);
 
 
        


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function huffman_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',false);
global image;
global filename;
global pathname;
a=image;
I=rgb2gray(a);
[m,n]=size(I);
Totalcount=m*n;
cnt=1;
sigma=0;
for i=0:255
k=I==i;
count(cnt)=sum(k(:));
pro(cnt)=count(cnt)/Totalcount;
sigma=sigma+pro(cnt);
cumpro(cnt)=sigma;
cnt=cnt+1;
end
symbols = [0:255];
dict = huffmandict(symbols,pro);
vec_size = 1;
for p = 1:m
for q = 1:n
newvec(vec_size) = I(p,q);
vec_size = vec_size+1;
end
end
hcode = huffmanenco(newvec,dict);
dhsig1 = huffmandeco(hcode,dict);
dhsig = uint8(dhsig1);
dec_row=sqrt(length(dhsig));
arr_row = 1;
arr_col = 1;
vec_si = 1;
for x = 1:m
for y = 1:n
back(x,y)=dhsig(vec_si);
arr_col = arr_col+1;
vec_si = vec_si + 1;
end
arr_row = arr_row+1;
end
[deco, map] = gray2ind(back,256);
RGB = ind2rgb(deco,map);
axes(handles.axes2);
imshow(RGB);
infoEski=imfinfo([pathname,filename]);
imwrite(RGB,filename,'jpg');
    infoYeni=imfinfo(filename);
    delete(filename);
    set(handles.edit3,'String',"Eski dosya boyutu:"+infoEski.FileSize+" Yeni dosya boyutu:"+infoYeni.FileSize);


% --------------------------------------------------------------------
function dct_Callback(hObject, eventdata, handles)
set(handles.slider1,'Visible',true);
set(handles.slider1,'max',100);
global filename;
global pathname;
global image;
A=image;
Q=[16 11 10 16 24 40 51 61; 12 12 14 19 26 58 60 55; 14 13 16 24 40 57 69 56; 14 17 22 29 51 87 80 62; 18 22 37 56 68 109 103 77; 24 35 55 64 81 104 113 92; 49 64 78 87 103 121 120 101; 72 92 95 98 112 100 103 99 ];
X=mat2cell(A,8*ones(1,32),8*ones(1,32));
for i=1:32
    for j=1:32
       X{i,j}=dct2(X{i,j});
    end
end
QF=get(handles.slider1,'Value');
if QF<50
    S=5000/QF;
else
    S=200-2*QF;
end
T=floor((S*Q)+50)/100;
T(T==0)=1;
for i=1:32
    for j=1:32
        X{i,j}=round(X{i,j}./T);
    end
end
for i=1:32
    for j=1:32
        X{i,j}=X{i,j}.*T;
    end
end
for i=1:32
    for j=1:32
       X{i,j}=idct2(X{i,j});
    end
end
CI=cell2mat(X);
CI=uint8(CI);
axes(handles.axes2);
imshow(CI);
infoEski=imfinfo([pathname,filename]);
imwrite(CI,filename,'jpg');
    infoYeni=imfinfo(filename);
    delete(filename);
    set(handles.edit3,'String',"Eski dosya boyutu:"+infoEski.FileSize+" Yeni dosya boyutu:"+infoYeni.FileSize);


 
