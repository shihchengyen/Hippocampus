function placebyview()
  
%     vms = vmsv('auto','NumShuffles',1000,'FiltLowOcc',0);
    vms = vmsv('auto');
    
%     vms = load('vmsv_filt.mat');
%     vms = vms.vms;
    
    vms = vms.data;
    map_used = vms.maps_raw;
%     vmp = vmpc('auto');
%     vmp = vmp.data;
    um = umaze('auto');
    um = um.data;
    st = load('spiketrain.mat');
    st = st.timestamps;
    st = st./1000;
    
    figure('name',vms.origin{1},'units','normalized','outerposition',[0.1 0.1 0.8 0.8])
    clf;
    
    binDepths = vms.binDepths;
    view_bin_count = 5122;
    index_tracker = 1:view_bin_count;

    % Restructure bins from linear to separate grids
    lin_grid_reference = cell(size(binDepths,1),1);

    for jj = 1:size(binDepths,1) % for each grid
        % Initialise empty matrices
        gridded_index = nan(binDepths(jj,1),binDepths(jj,2));

        % Assign linear bin to grid bin
        for mm = 1:binDepths(jj,1)*binDepths(jj,2) % For every point in linear map
            if mod(mm,binDepths(jj,2)) == 0
                y = binDepths(jj,2);
            else
                y = mod(mm,binDepths(jj,2));
            end
            x = ceil(mm/binDepths(jj,2));
            indbins_lin = mm + sum(binDepths(1:jj-1,1).*binDepths(1:jj-1,2));
            % Assign
            gridded_index(x,y) = index_tracker(indbins_lin);

        end
        % Collect output
        lin_grid_reference{jj} = gridded_index;
    end    
    
    dummy_pillars = lin_grid_reference;
    for jj = 6:9 % setting pillars to NaN
        dummy_pillars{jj} = NaN(size(dummy_pillars{jj}));
    end
   

    disp('calculating pair-wise sessionTime');
    
    temp_view = vms.sessionTime_generated(:,1:2);
    temp_view = [temp_view zeros(size(temp_view,1),1) ones(size(temp_view,1),1)];
    temp_place = um.sessionTime(:,1:2);
    temp_place = [temp_place(:,1) zeros(size(temp_place,1),1) temp_place(:,2) 2.*ones(size(temp_place,1),1)];
    temp_comb = [temp_view; temp_place];
    temp_comb = sortrows(temp_comb, 1);
    
    for row = 2:size(temp_comb,1)
        if temp_comb(row,4) == 1
            temp_comb(row,3) = temp_comb(row-1,3);
        else
            temp_comb(row,2) = temp_comb(row-1,2);
        end
    end
    
    spikings = histcounts(st,temp_comb(:,1)');
    spikings = [spikings 0]';
    
    temp_comb(1:end-1,4) = temp_comb(2:end,1) - temp_comb(1:end-1,1);
    temp_comb = [temp_comb spikings];
    occ_duration = zeros(view_bin_count, um.gridSteps*um.gridSteps);
    total_spikes = occ_duration;
    for row = 1:size(occ_duration,1)
        subset = temp_comb(temp_comb(:,2)==row,3:5);
        if ~isempty(subset)
            unique_place = unique(subset(:,1));
            unique_place(isnan(unique_place)) = [];
            unique_place(unique_place==0) = [];
            for place = 1:length(unique_place)
                occ_duration(row,unique_place(place)) = sum(subset(subset(:,1)==unique_place(place), 2));
                total_spikes(row,unique_place(place)) = sum(subset(subset(:,1)==unique_place(place), 3));
            end
        end
    end
    if vms.Args.FiltLowOcc
        occ_duration(sum(occ_duration>0,2)<=vms.MinTrials,:) = NaN;
    end
    occ_duration(occ_duration==0) = NaN;
    fr = total_spikes./occ_duration;
    
    h0 = axes('Position',[0.025 0.2 0.4 0.6],'Tag','h0_tag');

        c0 = colormap(h0, jet);
        colorbar(h0,'Position',[0.45 0.2 0.025 0.6]);
        for i = 3:9
            hold on;
            plotspatialview(i,map_used{i},lin_grid_reference{i},gca,0);
            hold off;
        end
        box off;

    set(h0,'UserData','rotate');
    cursorMode = rotate3d(gcf);
    set(cursorMode,'enable','on');
    hManager = uigetmodemanager(gcf);
    [hManager.WindowListenerHandles.Enabled] = deal(false);
   
    
    h2 = axes('Position',[0.525 0.2 0.40 0.6],'Tag','h2_tag');
    
            map = lin_grid_reference{3};
            map = rot90(map);

            surfx = repmat((0:40)',1,41);
            surfy = repmat(0:40,41,1);
            surfz = -1.*ones(41);            

            c2 = colormap(h2, gray);
            colorbar(h2,'Position',[0.95 0.2 0.025 0.6]);
            hold on;
            surf(h2, surfx,surfy,surfz,map);
            hold off;
            plotspatialview(6,0,0,gca,1);
            plotspatialview(7,0,0,gca,1);
            plotspatialview(8,0,0,gca,1);
            plotspatialview(9,0,0,gca,1);    
    
    linkprop([h0 h2], {'View', 'XLim', 'YLim', 'ZLim',});  %linkprop([h0 h2], / {'View', 'XLim', 'YLim', 'ZLim',});
    set(gcf,'WindowButtonMotionFcn',{@hovercallback, gcf, h0, h2});
    set(gcf,'WindowScrollWheelFcn',{@scrollcallback,gcf,h0,h2,fr});
    
end
    
function hovercallback(~,~,gcf,h0,h2)

    linkprop([h0 h2], {'View', 'XLim', 'YLim', 'ZLim',}); 

end

function scrollcallback(~,~,gcf,h0,h2,fr)
    
    disp('registered');
    if 1==1
        disp('changing cursor mode');
        if strcmpi(get(h0, 'UserData'), 'rotate')
            cursorMode = datacursormode(gcf);
            set(cursorMode, 'enable', 'on');
            set(cursorMode,'UpdateFcn',{@customize_label, cursorMode});
            set(h0, 'UserData', 'datapoint');
            hManager = uigetmodemanager(gcf);
            [hManager.WindowListenerHandles.Enabled] = deal(false); 
            set(gcf,'WindowButtonUpFcn',{@selectcallback, gcf, h0, h2, fr}); 
            linkprop([h0 h2], {'View', 'XLim', 'YLim', 'ZLim',}); 
        else
            cursorMode = rotate3d(gcf);
            set(cursorMode, 'enable', 'on');
            set(h0, 'UserData', 'rotate');
            linkprop([h0 h2], {'View', 'XLim', 'YLim', 'ZLim',}); 
        end   
        
        hManager = uigetmodemanager(gcf);
        [hManager.WindowListenerHandles.Enabled] = deal(false); 
        set(gcf,'WindowScrollWheelFcn',{@scrollcallback,gcf,h0, h2,fr});
        set(gcf,'WindowButtonMotionFcn',{@hovercallback, gcf, h0, h2});
    end

end

function [text] = customize_label(~, event_obj, dcm_obj)

        info_struct = getCursorInfo(dcm_obj);
        pos = info_struct.Position;
        x_ind = get(info_struct.Target,'XData')==pos(1);
        y_ind = get(info_struct.Target,'YData')==pos(2);
        z_ind = get(info_struct.Target,'ZData')==pos(3);
        loc = x_ind & y_ind & z_ind;

        fr = get(info_struct.Target,'CData');
        if ~(sum(loc(size(loc,1),:)) == 1 || sum(loc(:,size(loc,2))) == 1)
            [a, b] = ind2sub(size(loc),find(loc==1));
            firing_rate = fr(a,b);   
            text = ['firing rate: ' num2str(firing_rate)];
        end        

end

function selectcallback(h,ed,gcf,h0,h2,fr)

    linkprop([h0 h2], {'View'});
    current_axes = gca;
    
    if strcmpi(get(h0, 'UserData'), 'datapoint') && strcmpi(current_axes.Tag, 'h0_tag')
        
        dcm_obj = datacursormode(gcf);
        info_struct = getCursorInfo(dcm_obj);
        pos = info_struct.Position;
        x_ind = get(info_struct.Target,'XData')==pos(1);
        y_ind = get(info_struct.Target,'YData')==pos(2);
        z_ind = get(info_struct.Target,'ZData')==pos(3);
        loc = x_ind & y_ind & z_ind;
        c_ind = get(info_struct.Target,'CData');
        [a, b] = ind2sub(size(loc),find(loc==1));
        disp('debug');
        c_ind(a,b)
        
        ud = get(info_struct.Target,'UserData');
        if ~(sum(loc(size(loc,1),:)) == 1 || sum(loc(:,size(loc,2))) == 1)
            
            [a, b] = ind2sub(size(loc),find(loc==1));
            bin_number = ud(a,b);
         
            map = reshape(fr(bin_number,:),40,40);
            map = rot90(map);

            surfx = repmat((0:40)',1,41);
            surfy = repmat(0:40,41,1);
            surfz = -1.*ones(41);            
            
            axes(h2);
            cla(h2);

            surf(h2, surfx,surfy,surfz,map);
            plotspatialview(6,0,0,gca,1);
            plotspatialview(7,0,0,gca,1);
            plotspatialview(8,0,0,gca,1);
            plotspatialview(9,0,0,gca,1);

                c_temp = copper(64);
                c_temp = flipud(c_temp);
                c_temp(1,:) = [0.1425 0.7077 0.95];
                c2 = colormap(h2, c_temp);
                hold on;
                colorbar(h2,'Position',[0.95 0.2 0.025 0.6]);
                hold off;               
            
            if nanmax(nanmax(map)) > 0
                set(h2,'CLim',[0 nanmax(nanmax(map))]);
            else         
                set(h2,'CLim',[100 200]);
            end
             
            box off;
            linkprop([h0 h2], {'View', 'XLim', 'YLim', 'ZLim',});                   

        end
    end
        
       hManager = uigetmodemanager(gcf);
       [hManager.WindowListenerHandles.Enabled] = deal(false);
       set(gcf,'WindowButtonUpFcn',{@selectcallback, gcf, h0, h2, fr});
       set(gcf,'WindowScrollWheelFcn',{@scrollcallback,gcf, h0, h2,fr});
       set(gcf,'WindowButtonMotionFcn',{@hovercallback, gcf, h0, h2});
end





