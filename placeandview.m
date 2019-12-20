function placeandview()

%     vms = vmsv('auto','NumShuffles',1000,'FiltLowOcc',0);
    vms = vmsv('auto');
    
    vms = vms.data;
    um = umaze('auto');
    um = um.data;
    st = load('spiketrain.mat');
    st = st.timestamps;
    st = st./1000;
    
    figure('name',vms.origin{1},'units','normalized','outerposition',[0.1 0.1 0.8 0.8])
%     figure(); % for older vms files
%     a = load('bindepths.mat'); % for older vms files
%     vms.binDepths = a.binDepths; % for older vms files
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
    
    
    % sic modification starts here
    
                    comb.lin_grid_reference = lin_grid_reference;
                    comb.total_spikes = total_spikes;
                    comb.occ_duration = occ_duration;
                    comb.temp_comb = temp_comb;
        
                    for i = 1:9
                        temp = comb.lin_grid_reference{i};
                        temp = temp(:,:,1);
                        comb.lin_grid_reference{i} = temp;
                    end

                    truth = comb.total_spikes ./ comb.occ_duration;
                    truth1 = comb.total_spikes;
                    truth2 = comb.occ_duration;
                    
                    comb.total_spikes(isnan(comb.total_spikes)) = 0;
                    comb.occ_duration(isnan(comb.occ_duration)) = 0;

                    multi_dim = cell(2,9);
                    for i = 3:9

                        smallest_grid = min(min(comb.lin_grid_reference{i}));
                        largest_grid = max(max(comb.lin_grid_reference{i}));
                        extracted_views_spikes = comb.total_spikes(smallest_grid:largest_grid,:);
                        extracted_views_dur = comb.occ_duration(smallest_grid:largest_grid,:);  

                        lgr = comb.lin_grid_reference{i};
                        if lgr(2,1) - lgr(1,1) == 1
                            reshaped_spikes = reshape(extracted_views_spikes,size(lgr,1),size(lgr,2),40,40);
                            reshaped_dur = reshape(extracted_views_dur,size(lgr,1),size(lgr,2),40,40);
                        elseif lgr(1,2) - lgr(1,1) == 1
                            reshaped_spikes = reshape(extracted_views_spikes,size(lgr,2),size(lgr,1),40,40);
                            reshaped_dur = reshape(extracted_views_dur,size(lgr,2),size(lgr,1),40,40);
                        end

                        multi_dim{1,i} = reshaped_spikes; %permute(reshaped_spikes,[2 1 3 4]);
                        multi_dim{2,i} = reshaped_dur; %permute(reshaped_dur,[2 1 3 4]);

                    end

                    multi_rate = cell(1,9);
                    for i = 3:9
                        for j = 1:2
                            multi_dim{j,i} = convn(multi_dim{j,i},ones(3,3,1,1),'same');
                        end
                        multi_rate{1,i} = multi_dim{1,i}./multi_dim{2,i};
                    end

                    recon_fc = NaN(5122,1600);
                    recon_fd = NaN(5122,1600);
                    recon_fc(1:2,:) = total_spikes(1:2,:);
                    recon_fd(1:2,:) = occ_duration(1:2,:);
                    for i = 3:9
                        smallest_grid = min(min(comb.lin_grid_reference{i}));
                        largest_grid = max(max(comb.lin_grid_reference{i}));
                        recon_fc(smallest_grid:largest_grid,:) = reshape(multi_dim{1,i},largest_grid-smallest_grid+1,1600);
                        recon_fd(smallest_grid:largest_grid,:) = reshape(multi_dim{2,i},largest_grid-smallest_grid+1,1600);
                    end
        

    fr = recon_fc ./ recon_fd;
    fc = recon_fc;
                    
    % sic modification ends here

    
    h0 = axes('Position',[0.025 0.2 0.4 0.6],'Tag','h0_tag');

        c0 = colormap(h0, jet);
        colorbar(h0,'Position',[0.45 0.2 0.025 0.6]);
        for i = 3:9
            hold on;
            plotspatialview2(i,zeros(size(lin_grid_reference{i})),lin_grid_reference{i},gca,0,1);
            hold off;
        end
        hold on;
        plotspatialview2(0,0,0,gca,1,0);
        hold off;
        box off;    
    
    set(h0,'UserData','rotate');
    cursorMode = rotate3d(gcf);
    set(cursorMode,'enable','on');
    hManager = uigetmodemanager(gcf);
    [hManager.WindowListenerHandles.Enabled] = deal(false);        
        
    h2 = axes('Position',[0.525 0.2 0.40 0.6],'Tag','h2_tag');
    
            map = lin_grid_reference{3};
            map = map(:,:,1);
            map = nan(size(map));
            map = rot90(map);

            surfx = repmat((0:40)',1,41);
            surfy = repmat(0:40,41,1);
            surfz = -1.*ones(41);            

            c2 = colormap(h2, gray);
            colorbar(h2,'Position',[0.95 0.2 0.025 0.6]);
            hold on;
            surf(h2, surfx,surfy,surfz,map,'UserData',rot90(reshape(1:1600,40,40)));
            hold off;
            plotspatialview2(6,0,0,gca,1,0);
            plotspatialview2(7,0,0,gca,1,0);
            plotspatialview2(8,0,0,gca,1,0);
            plotspatialview2(9,0,0,gca,1,0); 
            plotspatialview2(0,0,0,gca,1,0);
            

%     linkprop([h0 h2], {'View', 'XLim', 'YLim', 'ZLim',});
    set(gcf,'WindowScrollWheelFcn',{@scrollcallback,gcf,h0,h2,fr,fc,lin_grid_reference, binDepths});    

end

function scrollcallback(~,~,gcf,h0,h2,fr,fc,lin_grid_reference, binDepths)
    
    disp('registered');

    if 1==1
        disp('changing cursor mode');
        if strcmpi(get(h0, 'UserData'), 'rotate')
            cursorMode = datacursormode(gcf);
            set(cursorMode, 'enable', 'on');
%             set(cursorMode,'UpdateFcn',{@customize_label, cursorMode});
            set(h0, 'UserData', 'datapoint');
            hManager = uigetmodemanager(gcf);
            [hManager.WindowListenerHandles.Enabled] = deal(false); 
            set(gcf,'WindowButtonUpFcn',{@selectcallback, gcf, h0, h2, fr, fc, lin_grid_reference, binDepths}); 
%             linkprop([h0 h2], {'View', 'XLim', 'YLim', 'ZLim',}); 
        else
            cursorMode = rotate3d(gcf);
            set(cursorMode, 'enable', 'on');
            set(h0, 'UserData', 'rotate');
%             linkprop([h0 h2], {'View', 'XLim', 'YLim', 'ZLim',}); 
        end   
        
        hManager = uigetmodemanager(gcf);
        [hManager.WindowListenerHandles.Enabled] = deal(false); 
        set(gcf,'WindowScrollWheelFcn',{@scrollcallback,gcf,h0, h2,fr,fc,lin_grid_reference, binDepths});
    end

end




function selectcallback(h,ed,gcf,h0,h2,fr,fc,lin_grid_reference, binDepths)

%     linkprop([h0 h2], {'View'});
    current_axes = gca;
    
    if strcmpi(get(h0, 'UserData'), 'datapoint') && strcmpi(current_axes.Tag, 'h0_tag')
              
        dcm_obj = datacursormode(gcf);
        info_struct = getCursorInfo(dcm_obj);
        
        check_if_compass_selected = get(info_struct.Target,'UserData');
        pass = 1;
        if (size(check_if_compass_selected,1) == 1) && (size(check_if_compass_selected,2) == 2)
            pass = 0;
        end
            
        if pass == 1
            pos = info_struct.Position;
            x_ind = get(info_struct.Target,'XData')==pos(1);
            y_ind = get(info_struct.Target,'YData')==pos(2);
            z_ind = get(info_struct.Target,'ZData')==pos(3);
            loc = x_ind & y_ind & z_ind;
            c_ind = get(info_struct.Target,'CData');
            [a, b] = ind2sub(size(loc),find(loc==1));

            ud = get(info_struct.Target,'UserData');
            if ~(sum(loc(size(loc,1),:)) == 1 || sum(loc(:,size(loc,2))) == 1)

                [a, b] = ind2sub(size(loc),find(loc==1));
                bin_number = ud(a,b);
                disp(bin_number);

                map = reshape(fr(bin_number,:),40,40);
                map = rot90(map);
                fc_map = reshape(fc(bin_number,:),40,40);
                fc_map = rot90(fc_map);            

                surfx = repmat((0:40)',1,41);
                surfy = repmat(0:40,41,1);
                surfz = -1.*ones(41);            

                axes(h2);

                surf(h2, surfx,surfy,surfz,map,'UserData',rot90(reshape(1:1600,40,40)));
                set(h2, 'Tag', 'h2_tag');
                plotspatialview2(6,0,0,gca,1,0);
                plotspatialview2(7,0,0,gca,1,0);
                plotspatialview2(8,0,0,gca,1,0);
                plotspatialview2(9,0,0,gca,1,0);
                plotspatialview2(0,0,0,gca,1,0);            

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
            end
        end
        
    elseif strcmpi(get(h0, 'UserData'), 'datapoint') && strcmpi(current_axes.Tag, 'h2_tag')
        
        disp('placemap selected');
        dcm_obj = datacursormode(gcf);
        info_struct = getCursorInfo(dcm_obj);
        
        check_if_compass_selected = get(info_struct.Target,'UserData');
        pass = 1;
        if (size(check_if_compass_selected,1) == 1) && (size(check_if_compass_selected,2) == 2)
            pass = 0;
        end
            
        if pass == 1
            pos = info_struct.Position;
            x_ind = get(info_struct.Target,'XData')==pos(1);
            y_ind = get(info_struct.Target,'YData')==pos(2);
            z_ind = get(info_struct.Target,'ZData')==pos(3);
            loc = x_ind & y_ind & z_ind;
            c_ind = get(info_struct.Target,'CData');
            [a, b] = ind2sub(size(loc),find(loc==1));

            ud = get(info_struct.Target,'UserData');
            if ~(sum(loc(size(loc,1),:)) == 1 || sum(loc(:,size(loc,2))) == 1)

                [a, b] = ind2sub(size(loc),find(loc==1));
                bin_number = ud(a,b);
                disp(bin_number);
                
                    split_plot = cell(size(binDepths,1),1);

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
                            gridded_index(x,y) = fr(indbins_lin,bin_number);

                        end
                        % Collect output
                        split_plot{jj} = gridded_index;
                    end    
                
                        axes(h0);
                        cla;
                        for i = 3:9
                            hold on;
                            plotspatialview2(i,split_plot{i},lin_grid_reference{i},gca,0,1);
                            hold off;
                        end
                        hold on;
                        plotspatialview2(0,0,0,gca,1,0);
                        hold off;
                        box off;     
                        
                        c_temp = copper(64);
                        c_temp = flipud(c_temp);
                        c_temp(1,:) = [0.1425 0.7077 0.95];
                        c0 = colormap(h0, c_temp);
                        hold on;
                        colorbar(h0,'Position',[0.45 0.2 0.025 0.6]);
                        hold off;    
                
            end
        end

        
    end
        
       hManager = uigetmodemanager(gcf);
       [hManager.WindowListenerHandles.Enabled] = deal(false);
       set(gcf,'WindowButtonUpFcn',{@selectcallback, gcf, h0, h2, fr,fc,lin_grid_reference, binDepths});
       set(gcf,'WindowScrollWheelFcn',{@scrollcallback,gcf, h0, h2,fr, fc,lin_grid_reference, binDepths});
end




