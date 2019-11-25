function varargout = huyang1(varargin)
% HUYANG1 MATLAB code for huyang1.fig
%      HUYANG1, by itself, creates a new HUYANG1 or raises the existing
%      singleton*.
%
%      H = HUYANG1 returns the handle to a new HUYANG1 or the handle to
%      the existing singleton*.
%
%      HUYANG1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HUYANG1.M with the given input arguments.
%
%      HUYANG1('Property','Value',...) creates a new HUYANG1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before huyang1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to huyang1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help huyang1

% Last Modified by GUIDE v2.5 23-Nov-2019 17:47:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @huyang1_OpeningFcn, ...
                   'gui_OutputFcn',  @huyang1_OutputFcn, ...
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


% --- Executes just before huyang1 is made visible.
function huyang1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to huyang1 (see VARARGIN)

cd('session01/array01/channel019');
rh = rplhighpass('auto');
ts1 = rh.data.analogData;
cd('../../../');
cd('sessioneye/array01/channel019');
rh = rplhighpass('auto');
ts2 = rh.data.analogData;
ts = [ts1; ts2];
cd('../../../');
cd('mountains/channel019/output')
out1 = readmda('firings.mda');
tt=[]
windowl=20
maxcluster=max(out1(3,:))
handles.peak=cell(1,maxcluster)
for k=1:maxcluster
    tt=[]
    for i=1:size(out1,2)
        if out1(3,i)==k
            tt=[tt, out1(2,i)];
        end
    end
        for i=1:numel(tt)
            peakv=0;
            peaki=0;
            if ts(tt(i))<0
                for j=(tt(i):tt(i)+windowl)
                    if ts(j)<peakv
                        peaki=j;
                        peakv=ts(j);
                    end
                end
                handles.peak{k}(1,i)=peaki;
                handles.peak{k}(2,i)=peakv;
                
            else
                for j=(tt(i):tt(i)+windowl)
                    if ts(j)>peakv
                        peakv=ts(j);
                        peaki=j;
                    end
                end
                handles.peak{k}(1,i)=peaki;
                handles.peak{k}(2,i)=peakv;
                
            end
        end
end

interval1=[]
interval=handles.peak{1}(1,:);
interval=sort(interval);
for i=2:size(interval,2)
    inter=interval(1,i)-interval(1,i-1);
    intertime=inter/30000;
    interval1=[interval1,intertime];
end

handles.maxcluster=maxcluster
handles.index=1

handles.clustertext=['cluster ' num2str(handles.index) ' of ' num2str(handles.maxcluster)]
set(handles.clusternumber,'string',handles.clustertext)

handles.readytoinput=0
handles.click=0
handles.xcheck=[]
handles.ycheck=[]
handles.xnew=[]
handles.ynew=[]

axes(handles.histo);
histogram(handles.peak{1}(2,:))
axes(handles.dot);
scatter(handles.peak{1}(1,:),handles.peak{1}(2,:))
axes(handles.interval);
histogram(interval1)

% Choose default command line output for huyang1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes huyang1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = huyang1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
if handles.index-1>0
    handles.index=handles.index-1
    interval1=[]
    interval=handles.peak{handles.index}(1,:);
    interval=sort(interval);
    for i=2:size(interval,2)
        inter=interval(1,i)-interval(1,i-1);
        intertime=inter/30000;
        interval1=[interval1,intertime];
    end
    axes(handles.dot);
    scatter(handles.peak{handles.index}(1,:),handles.peak{handles.index}(2,:))
    axes(handles.histo);
    histogram(handles.peak{handles.index}(2,:))
    axes(handles.interval);
    histogram(interval1)

    handles.clustertext=['cluster ' num2str(handles.index) ' of ' num2str(handles.maxcluster)]
    set(handles.clusternumber,'string',handles.clustertext)
     
end

guidata(hObject, handles);
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)

% for j=1:handles.maxcluster
%     for i=2:size(handles.peak{j},2)
%         inter=handles.peak{j}(1,i)-handles.peak{j}(1,i-1);
%         intertime=inter/30000;
%         handles.interval1{j}=[handles.interval1{j},intertime];
%     end
% end
if handles.index+1<handles.maxcluster+1
    interval1=[]
    handles.index=handles.index+1
    interval=handles.peak{handles.index}(1,:);
    interval=sort(interval);
    for i=2:size(interval,2)
        inter=interval(1,i)-interval(1,i-1);
        intertime=inter/30000;
        interval1=[interval1,intertime];
    end
    axes(handles.dot);
    scatter(handles.peak{handles.index}(1,:),handles.peak{handles.index}(2,:))
    axes(handles.histo);
    histogram(handles.peak{handles.index}(2,:))
    axes(handles.interval);
    histogram(interval1)
    
    handles.clustertext=['cluster ' num2str(handles.index) ' of ' num2str(handles.maxcluster)]
    set(handles.clusternumber,'string',handles.clustertext)
end

guidata(hObject, handles);
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function enternumber_Callback(hObject, eventdata, handles)

number=str2double(get(handles.enternumber,'string'))
if number<=handles.maxcluster & number>0
    handles.index=number
    interval1=[]
    interval=handles.peak{handles.index}(1,:);
    interval=sort(interval);
    for i=2:size(interval,2)
        inter=interval(1,i)-interval(1,i-1);
        intertime=inter/30000;
        interval1=[interval1,intertime];
    end
    axes(handles.dot);
    scatter(handles.peak{handles.index}(1,:),handles.peak{handles.index}(2,:))
    axes(handles.histo);
    histogram(handles.peak{handles.index}(2,:))
    axes(handles.interval);
    histogram(interval1)
   
    handles.clustertext=['cluster ' num2str(handles.index) ' of ' num2str(handles.maxcluster)]
    set(handles.clusternumber,'string',handles.clustertext)

end
guidata(hObject, handles);
% hObject    handle to enternumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enternumber as text
%        str2double(get(hObject,'String')) returns contents of enternumber as a double


% --- Executes during object creation, after setting all properties.
function enternumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enternumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select.
function select_Callback(hObject, eventdata, handles)
handles.readytoinput=1
guidata(hObject, handles);
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
handles.xaxis=handles.peak{handles.index}(1,:)
handles.yaxis=handles.peak{handles.index}(2,:)
handles.xcheck=zeros([1 size(handles.xaxis,2)])
handles.ycheck=zeros([1 size(handles.yaxis,2)])
if handles.click==1
    out = get(gca,'CurrentPoint');
    handles.index1=handles.index
    for i=1:size(handles.xaxis,2)
        if handles.xaxis(i)<out(1,1)+500000 & handles.xaxis(i)>out(1,1)-500000
            if handles.yaxis(i)<out(1,2)+2 & handles.yaxis(i)>out(1,2)-2
               handles.xcheck(i)=1
               handles.ycheck(i)=1
             
%                hold on
%                axes(handles.dot);
%                scatter(handles.xaxis(i),handles.yaxis(i),'r')
%                hold off
            end
        end
    end
    for i = 1:size(handles.xcheck,2)
        if handles.xcheck(i)==1
            handles.xnew=[handles.xnew,handles.xaxis(i)]
            handles.ynew=[handles.ynew,handles.yaxis(i)]
            handles.xnew=unique(handles.xnew,'stable');
            handles.ynew=unique(handles.ynew,'stable');
        end
    end
   
end




guidata(hObject, handles);

% --- Executes on button press in finish.
function finish_Callback(hObject, eventdata, handles)
 handles.peak{handles.maxcluster+1}=[handles.xnew;handles.ynew]
 handles.maxcluster=handles.maxcluster+1
 
 hold on
 axes(handles.dot);
 scatter(handles.xnew,handles.ynew,'r')
 hold off
 
 handles.index1text=['from' num2str(handles.index1) ]
 set(handles.from,'string',handles.index1text)

 handles.xnew=[]
 handles.ynew=[]
 handles.xaxis=[]
 handles.yaxis=[]
 handles.xcheck=[]
 handles.ycheck=[]
 handles.clustertext=['cluster ' num2str(handles.index) ' of ' num2str(handles.maxcluster)]
 
 set(handles.clusternumber,'string',handles.clustertext)
 guidata(hObject, handles);
 

% hObject    handle to finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
if handles.readytoinput==1
    if handles.click==0
        handles.click=1
    else
        handles.click=0
        handles.readytoinput=0
    end
else 
    handles.click=0
end
guidata(hObject, handles);
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
