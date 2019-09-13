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

    rough = readmda('firings.curated.mda');
    time_series = readmda('dataset/raw_data.mda');

    unique_cells_here = unique(rough(3,:));
    handles.cell_count = length(unique_cells_here);
    handles.cell_list = NaN(1,100);
    handles.cell_list(1:length(unique_cells_here)) = unique_cells_here;
    handles.full_data = NaN(152, length(rough(1,:)));
    handles.full_data(1:2,:) = rough(2:3,:);
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
    
    handles.initial_large_clusters = NaN(handles.cell_count, length(handles.full_data(1,:)) + 2);
    
    fat_index = 1;
    for cell_index = 1:length(unique_cells_here)
        [~, col] = find(handles.full_data(2,:)==unique_cells_here(cell_index));
        if length(col) > 15000
            handles.initial_large_clusters(fat_index, 1) = unique_cells_here(cell_index);
            handles.initial_large_clusters(fat_index, 2) = length(col);
            handles.initial_large_clusters(fat_index, 5:4+length(col)) = col;
            fat_index = fat_index + 1;
        end
    end
 
    base = gcf();
    
    [a, b] = unix('find . -name "tmp" -maxdepth 1 | wc -l');
    disp(b);

    if str2double(b) == 0
        mkdir('tmp');
        cd('tmp');
        disp(handles.initial_large_clusters(1:fat_index-1,1:2));
        for large_clusters_index = 1:fat_index-1
            cluster_name = handles.initial_large_clusters(large_clusters_index, 1);
            temp_data = handles.full_data(3:152,handles.initial_large_clusters(large_clusters_index,5:4+handles.initial_large_clusters(large_clusters_index, 2)));
            handles.initial_large_clusters(large_clusters_index, 4) = max(max(temp_data));
            handles.initial_large_clusters(large_clusters_index, 3) = min(min(temp_data));
            image_count = 1;
            reading_anchor = 1;         
            while reading_anchor < handles.initial_large_clusters(large_clusters_index, 2)
                indices = handles.initial_large_clusters(large_clusters_index, 5:4+handles.initial_large_clusters(large_clusters_index, 2));
                to_use = indices(reading_anchor:min(reading_anchor+9999, length(indices)));
                plotting_data = handles.full_data(3:152,to_use);

                disp(size(plotting_data));
                disp(reading_anchor);

                fig2 = figure();
                ax1 = axes('Parent', fig2); 
                plot(ax1, [1:150], plotting_data, 'Color', [0 0 0], 'LineWidth', 0.25);
                xlim([1 150]);
                ylim(handles.initial_large_clusters(large_clusters_index, 3:4));            
                set(ax1, 'Color', 'None');

                    outerpos = ax1.OuterPosition;
                    ti = ax1.TightInset; 
                    left = outerpos(1) + ti(1);
                    bottom = outerpos(2) + ti(2);
                    ax_width = outerpos(3) - ti(1) - ti(3);
                    ax_height = outerpos(4) - ti(2) - ti(4);
                    ax1.Position = [left bottom ax_width ax_height];

                set(ax1,'XTickLabel',[], 'YTickLabel', []);
                set(ax1, 'Color', 'None');                    
                
                set(ax1, 'XLimSpec', 'Tight');
                set(ax1, 'YLimSpec', 'Tight');            
                
                file_name1 = char(strcat(string(cluster_name),'_', string(image_count), '.png'));
                export_fig(file_name1, '-transparent', '-png');
                clf(fig2, 'reset');
                close(fig2);

                reading_anchor = reading_anchor + 10000;
                image_count = image_count + 1;
            end
        end    
        to_save = handles.initial_large_clusters;
        save('table.mat', 'to_save');
        cd('..');
    else
        disp('no need to regenerate');
        handles.initial_large_clusters = load('tmp/table.mat');
        handles.initial_large_clusters = handles.initial_large_clusters.to_save;
    end

    handles.view_index = 1;
    handles.number_selected = 0;
    disp(handles.initial_large_clusters(:,1:4));
    figure(base);
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
varargout{1} = handles.output;


function [handles] = update_plot(hObject, eventdata, handles)

    disp('called');
    handles.target = handles.cell_list(handles.view_index);
    link_subset_to_full = NaN(1,size(handles.full_data,2));
    possible_waves = NaN(150,size(handles.full_data,2));
    temp_index = 1;
    for col = 1:size(handles.full_data,2)
        if handles.full_data(2,col) == handles.target
            possible_waves(:,temp_index) = handles.full_data(3:152,col);
            link_subset_to_full(1,temp_index) = col;
            temp_index = temp_index + 1;
        end
    end
    link_subset_to_full = link_subset_to_full(1:temp_index-1);
    possible_waves = possible_waves(:,1:temp_index-1);
    cla(handles.axes1, 'reset');
    set(handles.cluster_id, 'String', 'Current Cluster ID');
    set(handles.selection_count, 'String', 'Count');
    
    handles.link_subset_to_full = link_subset_to_full;
    handles.possible_waves = possible_waves;
    handles.hit = zeros(1,size(handles.possible_waves,2));
    disp(handles.view_index);
    disp(handles.target);
    
    command = ['find ./tmp/ -name "', num2str(handles.target), '_*" | wc -l'];
    [~,out] = unix(command);
    disp('pngs found');
    disp(out);
    if str2double(out) > 0
        command = ['find ./tmp/ -name "', num2str(handles.target), '_*"'];
        [~,list] = unix(command);
        split = strsplit(list, './tmp//');
        split = split(2:length(split));
        for i = 1:length(split)
            title = split{i};
            title = title(1:length(title)-1);
            if i == 1
                xlim([1 150]);
                loc = find(handles.initial_large_clusters(:,1)==handles.target);
                disp(loc);
                disp(handles.initial_large_clusters(loc, 3:4));
                ylim(handles.initial_large_clusters(loc, 3:4));              
            end
            cd('tmp');
            [png1, png2, png3] = imread(title);    
            image(handles.axes2, png1, 'AlphaData', png3);
            set(handles.axes2,'XTickLabel',[], 'YTickLabel', []);            
            set(handles.axes2, 'XLimSpec', 'Tight');
            set(handles.axes2, 'YLimSpec', 'Tight');            
            if i == 1
                hold(handles.axes2, 'on');
            end
            cd('..');
        end
        hold(handles.axes2, 'off');
    else
        disp(size(handles.possible_waves));
        cla(handles.axes2, 'reset');
        handles.line_handles = plot(handles.axes1, 1:150, handles.possible_waves', 'Color', [0 0 0], 'LineWidth', 0.25);
    end
    
    found_for_range = 0;
    for i = 1:length(handles.initial_large_clusters(:,1))
        if handles.initial_large_clusters(i,1) == handles.target
            found_for_range = i;
            break;
        end
    end
    
    if found_for_range ~= 0
        ylim(handles.axes1, [handles.initial_large_clusters(found_for_range,3), handles.initial_large_clusters(found_for_range,4)]);
    end
    xlim(handles.axes1, [1 150]);  
    set(handles.axes1, 'Color', 'None');
    
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
    set(handles.selection_count, 'String', strcat(num2str(sum(handles.hit)), '/', num2str(size(handles.possible_waves, 2))));

guidata(hObject, handles);


% --- Executes on button press in pb7.
function [handles] = pb7_Callback(hObject, eventdata, handles)

    set(handles.pb7,'visible','off');
    set(handles.pb6,'visible','off');
    set(handles.pb5,'visible','off');
    set(handles.pb4,'visible','off');
    set(handles.pb3,'visible','off');
    set(handles.pb2,'visible','off');
    set(handles.pb1,'visible','off');
    
    
    [x,y] = ginput(2);
   
    if sum(handles.hit) > 0
        delete(handles.extra_lines);
    end
    
    for col = 1:size(handles.possible_waves,2)
        for points = 1:size(handles.possible_waves,1)
            if (min(y) < handles.possible_waves(points,col)) && (handles.possible_waves(points,col) < max(y))
                if (min(x) < points) && (points < max(x))
                    handles.hit(1,col) = 1;
                    break;
                end
            end
        end
    end
    [~, coli] = find(handles.hit);
    redraw = handles.possible_waves(:,coli);      
    
    hold(handles.axes1, 'on');
    if sum(handles.hit) > 0
        handles.extra_lines = plot(handles.axes1, 1:150, redraw', 'Color', [1 0 0], 'LineWidth', 0.25);
    end
    hold(handles.axes1, 'off');
    set(handles.selection_count, 'String', strcat(num2str(sum(handles.hit)), '/', num2str(size(handles.possible_waves, 2))));
    
    set(handles.pb7,'visible','on');
    set(handles.pb6,'visible','on');
    set(handles.pb5,'visible','on');
    set(handles.pb4,'visible','on');
    set(handles.pb3,'visible','on');
    set(handles.pb2,'visible','on');
    set(handles.pb1,'visible','on');
    


guidata(hObject, handles);
















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
