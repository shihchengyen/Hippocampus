function varargout = refinement(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @refinement_OpeningFcn, ...
                   'gui_OutputFcn',  @refinement_OutputFcn, ...
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


% --- Executes just before refinement is made visible.
function refinement_OpeningFcn(hObject, eventdata, handles, varargin)
    
    set(gcf, 'WindowButtonMotionFcn', '')
    set(handles.axes4, 'XTick', []);
    set(handles.axes4, 'YTick', []);
    handles.start_patch = 0;
    handles.corner_count = 0;
    set(handles.axes4, 'Color', 'None');
    
    handles.waiting_for_inputs = 0;
    rough = readmda('firings.curated.mda');
    if length(rough) == 0
        disp('no preliminary cells marked out for refinement, exiting.');
        figure1_CloseRequestFcn(hObject, eventdata, handles);
        return;
    end
    time_series = readmda('dataset/raw_data.mda');

    unique_cells_here = unique(rough(3,:));
    handles.cell_count = length(unique_cells_here);
    handles.cell_list = NaN(1,100);
    handles.cell_list(1:length(unique_cells_here)) = unique_cells_here;
    handles.full_data = NaN(153, length(rough(1,:)));
    handles.full_data(1:2,:) = rough(2:3,:); % first row sample number for spike occurrence, second row current cluster number
    handles.full_data(153,:) = handles.full_data(2,:); % last section to keep track of original grouping, for drawing over with dotted
    for col = 1:size(handles.full_data, 2)
        if handles.full_data(1,col)-74 < 1
            short = time_series(1:handles.full_data(1,col)+75);
            padded = padarray(short, [0, 150-length(short)], 'pre');
            handles.full_data(3:152,col) = padded';
        elseif handles.full_data(1,col)+75 > length(time_series)
            short = time_series(handles.full_data(1,col)-74:length(time_series));
            padded = padarray(short, [0, 150-length(short)], 'post');
            handles.full_data(3:152,col) = padded';
        else
            handles.full_data(3:152,col) = time_series(handles.full_data(1,col)-74:handles.full_data(1,col)+75)';
        end
    end
    
    handles.view_index = 1;
    handles.number_selected = 0;
    handles.last_index = length(time_series);
    
    [handles] = update_plot(hObject, eventdata, handles);

% Choose default command line output for refinement
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes refinement wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = refinement_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
    varargout{1} = handles.output;
catch
    disp('done.');
end

function [handles] = update_plot(hObject, eventdata, handles)

    disp('called');
    
    handles.target = handles.cell_list(handles.view_index);
    link_subset_to_full = NaN(1,size(handles.full_data,2));
    possible_waves = NaN(150,size(handles.full_data,2));
    amp_with_indices = NaN(2,size(handles.full_data,2));

    temp_index = 1;
    for col = 1:size(handles.full_data,2)
        if handles.full_data(2,col) == handles.target
            possible_waves(:,temp_index) = handles.full_data(3:152,col);
            link_subset_to_full(1,temp_index) = col;
            amp_with_indices(1,temp_index) = handles.full_data(1,col);
            amp_with_indices(2,temp_index) = handles.full_data(77,col);
            temp_index = temp_index + 1;
        end
    end
    link_subset_to_full = link_subset_to_full(1:temp_index-1);
    possible_waves = possible_waves(:,1:temp_index-1);
       
    set(handles.cluster_id, 'String', 'Current Cluster ID');
    set(handles.selection_count, 'String', 'Count');
    
    handles.link_subset_to_full = link_subset_to_full;
    handles.possible_waves = possible_waves;
    handles.amp_with_indices = amp_with_indices;
    
    handles.hit = zeros(1,size(handles.possible_waves,2));
    disp(handles.view_index);
    disp(handles.target);

        disp(size(handles.possible_waves));
%         cla(handles.axes1, 'reset');
        
        
        vals = NaN(151, size(handles.possible_waves, 2));
        vals(1:150,:) = handles.possible_waves;
        ts = ([1:150 150])'*ones(1,length(handles.possible_waves(1,:)));
        
        
        xlim(handles.axes1, [1 150]);
        ylim(handles.axes1, [min(min(handles.possible_waves)) max(max(handles.possible_waves))]);
        disp('mass plot start');
        plot(handles.axes1, ts(:), vals(:), 'Color', [0 0 0], 'LineWidth', 0.25);
        disp('mass plot end');
        
                          set(handles.axes1, 'XTick', []);
            set(handles.axes1, 'XLimSpec', 'Tight');
            set(handles.axes1, 'YLimSpec', 'Tight');   


%         cla(handles.axes3, 'reset');
        
        plot(handles.axes3, handles.amp_with_indices(1,:), handles.amp_with_indices(2,:), 'LineStyle', 'none', 'Marker', 'x', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        xlim(handles.axes3, [1 handles.last_index]);
        ylim(handles.axes3, [min(handles.amp_with_indices(2,:))-5 max(handles.amp_with_indices(2,:))+5]);
            
            set(handles.axes3, 'XLimSpec', 'Tight');
            set(handles.axes3, 'YLimSpec', 'Tight');   
            space = ' ';
            title1 = ['cluster', space, num2str(handles.view_index), space, 'of', space, num2str(handles.cell_count)];
            title(handles.axes3, title1(:)', 'FontSize', 14);
  
            
    set(handles.cluster_id, 'String', strcat('Cluster ID: ',num2str(handles.target)));
    set(handles.selection_count, 'String', strcat('0/', num2str(size(handles.possible_waves, 2))));

    disp('ended');
    
guidata(hObject, handles);

function [handles] = refresh_after_changes(hObject, eventdata, handles)

    unique_cells_here = unique(handles.full_data(2,:));
    disp(unique_cells_here);
    for i = 1:length(unique_cells_here)
        if unique_cells_here(i) == 0
            unique_cells_here(i) = [];
            break;
        end
    end
    
    handles.cell_count = length(unique_cells_here);
    handles.cell_list = NaN(1,100);
    handles.cell_list(1:length(unique_cells_here)) = unique_cells_here;
    
    handles.view_index = 1;
    handles.number_selected = 0;
    [handles] = update_plot(hObject, eventdata, handles);


guidata(hObject, handles);

% --- Executes on button press in pb2.
function [handles] = pb2_Callback(hObject, eventdata, handles)

    if sum(handles.hit) > 0
        
        for i = 1:length(handles.hit)
            if handles.hit(i) == 1
                link_back = handles.link_subset_to_full(i);
                handles.full_data(2,link_back) = 0;
            end
        end
        
        [handles] = refresh_after_changes(hObject, eventdata, handles);
        
    end


guidata(hObject, handles);


% --- Executes on button press in pb1.
function [handles] = pb1_Callback(hObject, eventdata, handles)

    if sum(handles.hit) > 0
        
        for i = 1:length(handles.hit)
            if handles.hit(i) == 1
                link_back = handles.link_subset_to_full(i);
                check_allocation = 0;
                for j = 1:length(handles.cell_list)
                    if (handles.cell_list(j) >= 100*handles.target) && (handles.cell_list(j) < 100*(handles.target+1))
                        check_allocation = check_allocation + 1;
                    end
                end
                allocated = (100*handles.target) + check_allocation;
                handles.full_data(2,link_back) = allocated;
            end
        end
        
        [handles] = refresh_after_changes(hObject, eventdata, handles);
        
    end

guidata(hObject, handles);

% --- Executes on button press in pb4.
function pb4_Callback(hObject, eventdata, handles)

    date_done = date;
    folder_name = strcat('prepostcuration_', date_done);
    s = ' ';
    unix(['mkdir', s, folder_name]);
    cd(folder_name);
    unix('cp ../firings.curated.mda .');
    base = readmda('firings.curated.mda');
    base(3,:) = handles.full_data(2,:);
    writemda(base, 'firings.curated2.mda');
    cd('..');
    base1 = base(:,base(3,:)>0);
    export_mountain_cells(base1);

    
delete(hObject);


% --- Executes on button press in pb6.
function [handles] = pb6_Callback(hObject, eventdata, handles)

    handles.hit = zeros(size(handles.hit));
    delete(handles.extra_lines);
        delete(handles.extra_dots);
        delete(handles.extra_dots2);   
%         delete(handles.anchor_point);
    set(handles.selection_count, 'String', strcat(num2str(sum(handles.hit)), '/', num2str(size(handles.possible_waves, 2))));

guidata(hObject, handles);

function [handles] = inputs_ready(hObject, eventdata, handles)

    x = handles.corner_data(:,1);
    y = handles.corner_data(:,2);

    if sum(handles.hit) > 0
        delete(handles.extra_lines);
        delete(handles.extra_dots);
        delete(handles.extra_dots2);
%         delete(handles.anchor_point);
    end
    left = round(min(x));
    right = round(max(x));
    up = round(max(y));
    down = round(min(y));
    
    if strcmp(handles.selected_axes, 'axes1')

        for col = 1:size(handles.possible_waves,2) % optimized version here, minimize rectangle width for faster searching
            if sum(down<handles.possible_waves(left:right,col) & handles.possible_waves(left:right,col)<up) > 0
                handles.hit(1,col) = 1;
            end
        end

    else
        
        logic_array1 = (handles.amp_with_indices(1,:) < right) & (handles.amp_with_indices(1,:) > left);
        logic_array2 = (handles.amp_with_indices(2,:) < up) & (handles.amp_with_indices(2,:) > down);
        
        handles.hit(1,logic_array1 & logic_array2) = 1;
        
    end
        
        
    disp(sum(handles.hit));
    [~, coli] = find(handles.hit);
    redraw = handles.possible_waves(:,coli);      
    
    disp('plotting extra');
    hold(handles.axes1, 'on');
    
        vals = NaN(151, size(redraw, 2));
        vals(1:150,:) = redraw;
        ts = ([1:150 150])'*ones(1,length(redraw(1,:)));
        
        disp('mass plot start');
        handles.extra_lines = plot(handles.axes1, ts(:), vals(:), 'Color', [1 0 0], 'LineWidth', 0.25);
        disp('mass plot end');
    
    hold(handles.axes1, 'off');
    disp('plotted extra');
    
    hold(handles.axes3, 'on');
    rescatter = handles.amp_with_indices(:,coli);
    handles.extra_dots = plot(handles.axes3, rescatter(1,:), rescatter(2,:), 'LineStyle', 'none', 'Marker', 'x', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    handles.extra_dots2 = plot(handles.axes3, rescatter(1,:), rescatter(2,:), 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r');
    hold(handles.axes3, 'off');
    
    set(handles.selection_count, 'String', strcat(num2str(sum(handles.hit)), '/', num2str(size(handles.possible_waves, 2))));
    

    set(handles.pb7,'visible','on');
    set(handles.pb6,'visible','on');
    set(handles.pb5,'visible','on');
    set(handles.pb4,'visible','on');
    set(handles.pb3,'visible','on');
    set(handles.pb2,'visible','on');
    set(handles.pb1,'visible','on');    


guidata(hObject, handles);


% --- Executes on button press in pb7.
function [handles] = pb7_Callback(hObject, eventdata, handles)

    tic;
    handles.waiting_for_inputs = 1;
    handles.corner_count = 0;
    handles.corner_data = NaN(2,2);
    disp(toc);
    
    tic;
%     set(handles.pb7, 'String', 'please select');
    set(handles.pb7,'visible','off');
    set(handles.pb6,'visible','off');
    set(handles.pb5,'visible','off');
    set(handles.pb4,'visible','off');
    set(handles.pb3,'visible','off');
    set(handles.pb2,'visible','off');
    set(handles.pb1,'visible','off');
    disp(toc);
    
    tic;
guidata(hObject, handles);
disp(toc);















% --- Executes on button press in pb3.
function pb3_Callback(hObject, eventdata, handles)

    if handles.view_index > 1
        handles.view_index = handles.view_index - 1;
        [handles] = update_plot(hObject, eventdata, handles);
    end

guidata(hObject, handles);

% --- Executes on button press in pb5.
function pb5_Callback(hObject, eventdata, handles)

    if handles.view_index < handles.cell_count
        handles.view_index = handles.view_index + 1;
        [handles] = update_plot(hObject, eventdata, handles);
    end

guidata(hObject, handles);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function [handles] = figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    disp('click');
    if handles.waiting_for_inputs ~= 1
        disp('kicked');
        return;
    end
    disp('passed');

    fig1 = gcf;
%     disp(fig1.SelectionType);
    if strcmp(fig1.SelectionType, 'alt')
                handles.waiting_for_inputs = 0;
                handles.corner_count = 0;
                handles.selected_axes = 'none';
%                 delete(handles.anchor_point);
                set(handles.pb7,'visible','on');
                set(handles.pb6,'visible','on');
                set(handles.pb5,'visible','on');
                set(handles.pb4,'visible','on');
                set(handles.pb3,'visible','on');
                set(handles.pb2,'visible','on');
                set(handles.pb1,'visible','on');   
                guidata(hObject, handles);
                return;
    end
    
    
    set(handles.axes1, 'Tag', 'axes1');
    set(handles.axes3, 'Tag', 'axes3');
    check1 = gca;
    check1.Tag
    
    if strcmp(check1.Tag, 'axes1') || strcmp(check1.Tag, 'axes3')
        disp('selected axes');
        if handles.corner_count == 0
            handles.selected_axes = check1.Tag;
        elseif strcmp(handles.selected_axes, check1.Tag) ~= 1
            return;
        end
    else
        return;
    end

    

            disp('getting point');
            out = get(gca,'CurrentPoint');
            disp('point received');

            handles.corner_count = handles.corner_count + 1;
            handles.corner_data(handles.corner_count,:) = out(1,1:2);
%             hold(check1, 'on');
%             handles.anchor_point = plot(check1, handles.corner_data(1,1), handles.corner_data(1,2), 'Marker', 'p', 'MarkerFaceColor', 'g', 'MarkerSize', 12);
%             hold(check1, 'off');
            
            if handles.corner_count == 2

                [handles] = inputs_ready(hObject, eventdata, handles);
                handles.waiting_for_inputs = 0;
                handles.corner_count = 0;
                handles.selected_axes = 'none';

            end



guidata(hObject, handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

