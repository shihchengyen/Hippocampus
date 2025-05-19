function [maps_adsm, durs_adsm, rad_adsm, maps_bcsm, maps_dksm, durs_bcsm, durs_dksm, rad_adsm_grid] = smoothMaps(maps_raw, dur_raw, spk_raw, spikes_count, Args)
                
                % Adaptive smoothing function for place maps
                durs_raw = repmat(dur_raw',1,Args.NumShuffles+1);
                preset_to_zeros = durs_raw == 0;
                
                % Switch from linear maps to grid maps
                durs_raw_grid = cell2mat(lineartogrid(durs_raw,'place',[Args.GridSteps Args.GridSteps]));
                spikes_count_grid = cell2mat(lineartogrid(spikes_count','place',[Args.GridSteps Args.GridSteps]));
                preset_to_zeros_grid = logical(cell2mat(lineartogrid(preset_to_zeros,'place',[Args.GridSteps Args.GridSteps])));
                maps_raw_grid = cell2mat(lineartogrid(maps_raw','place',[Args.GridSteps Args.GridSteps]));
                
                unvis = ~(durs_raw_grid>0);
                % Boxcar smoothing
                maps_bcsm_grid = smooth(maps_raw_grid,5,unvis,'boxcar');
                durs_bcsm_grid = smooth(durs_raw_grid,5,unvis,'boxcar');
                % Disk smoothing
                maps_dksm_grid = smooth(maps_raw_grid,5,unvis,'disk');
                durs_dksm_grid = smooth(durs_raw_grid,5,unvis,'disk');
               
                % Set up adaptive smoothing parameters and output vars
                to_compute = 1:0.5:Args.GridSteps/2;
                possible = NaN(length(to_compute),2,Args.GridSteps,Args.GridSteps,Args.NumShuffles + 1);
                maps_adsm_grid = NaN(size(possible,3), size(possible,4), size(possible,5));
                maps_adsm_grid(preset_to_zeros_grid) = 0;
                durs_adsm_grid = NaN(size(possible,3), size(possible,4), size(possible,5));
                durs_adsm_grid(preset_to_zeros_grid) = 0;
                rad_adsm_grid = NaN(size(possible,3), size(possible,4), size(possible,5));
                rad_adsm_grid(preset_to_zeros_grid) = 0;
                
                wip = ones(Args.NumShuffles+1,1);
                % Adaptive smoothing
                for idx = 1:length(to_compute)

                    f=fspecial('disk',to_compute(idx));
                    f(f>=(max(max(f))/3))=1;
                    f(f~=1)=0;

                    possible(idx,1,:,:,:) = repmat(imfilter(durs_raw_grid(:,:,1), f, 'conv'), 1,1,Args.NumShuffles+1);   %./scaler;
                    possible(idx,2,:,:,find(wip)) = imfilter(spikes_count_grid(:,:,find(wip)), f, 'conv');   %./scaler;

                    alphaValue = Args.Alpha;
                    logic1 = squeeze(alphaValue./(possible(idx,1,:,:,:).*sqrt(possible(idx,2,:,:,:))) <= to_compute(idx));
                    slice1 = squeeze(possible(idx,1,:,:,:));
                    slice2 = squeeze(possible(idx,2,:,:,:));

                    maps_adsm_grid(logic1 & isnan(maps_adsm_grid)) = slice2(logic1 & isnan(maps_adsm_grid))./slice1(logic1 & isnan(maps_adsm_grid));
                    durs_adsm_grid(logic1 & isnan(durs_adsm_grid)) = slice1(logic1 & isnan(durs_adsm_grid));
                    rad_adsm_grid(logic1 & isnan(rad_adsm_grid)) = to_compute(idx);

%                     disp('smoothed with kernel size:');
%                     disp(to_compute(idx));
%                     disp('grids left');
%                     disp(sum(sum(sum(isnan(to_fill(:,:,:))))));

                    check = squeeze(sum(sum(isnan(maps_adsm_grid),2),1));
                    wip(check==0) = 0;

                end
                
                % Reshape from grid to linear maps
                maps_adsm_grid(preset_to_zeros_grid) = nan; % unvisited bins should be nan
                maps_adsm = gridtolinear({maps_adsm_grid},'place',[Args.GridSteps Args.GridSteps]);
                maps_adsm = maps_adsm';
                durs_adsm_grid(isnan(durs_adsm_grid) | preset_to_zeros_grid) = 0;
                durs_adsm = gridtolinear({durs_adsm_grid},'place',[Args.GridSteps Args.GridSteps]);
                durs_adsm = durs_adsm';
                rad_adsm_grid(preset_to_zeros_grid) = nan;
                rad_adsm = gridtolinear({rad_adsm_grid},'place',[Args.GridSteps Args.GridSteps]);
                rad_adsm = rad_adsm';
                maps_bcsm = gridtolinear({maps_bcsm_grid},'place',[Args.GridSteps Args.GridSteps]);
                maps_bcsm = maps_bcsm';
                maps_dksm = gridtolinear({maps_dksm_grid},'place',[Args.GridSteps Args.GridSteps]);
                maps_dksm = maps_dksm';
                durs_bcsm_grid(isnan(durs_bcsm_grid) | preset_to_zeros_grid) = 0;
                durs_bcsm = gridtolinear({durs_bcsm_grid},'place',[Args.GridSteps Args.GridSteps]);
                durs_bcsm = durs_bcsm';
                durs_dksm_grid(isnan(durs_dksm_grid) | preset_to_zeros_grid) = 0;
                durs_dksm = gridtolinear({durs_dksm_grid},'place',[Args.GridSteps Args.GridSteps]);
                durs_dksm = durs_dksm';
end
