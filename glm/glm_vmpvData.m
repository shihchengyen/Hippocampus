function vmpvData = glm_vmpvData(tbin_size)
%
%   Run in directory of spiketrain.mat, will return a structure
%   
%   Used to bin vmpv data by time, to be used for glm fitting.
%   The smaller tbin_size is, the more accurate the place/view of allocated
%   spikes.
%
%   vmpvData.bin_stc contains the following:
%   column 1 - start of time bin
%   column 2 - place bin no.
%   column 3 - head direction bin no.
%   column 4 - view bin no.
%   column 5 - spikes occurred
%

pv = vmpv('auto');

spiketrain = load('spiketrain.mat');
spiketimes = spiketrain.timestamps/1000; % now in seconds
stc = pv.data.sessionTimeC;

% using vel filter only for now
ThresVel = 1;
conditions = ones(size(stc,1),1);
conditions = conditions & get(pv,'SpeedLimit',ThresVel);

% don't touch for now
UseMinObs = false;
if UseMinObs
    place_bins_sieved = pv.data.place_good_bins;
    view_bins_sieved = pv.data.view_good_bins;
    conditions = conditions & (pv.data.pv_good_rows);
else
    place_bins_sieved = 1:(40 * 40);
    view_bins_sieved = 1:5122;
end

% Construct new stc, with each row representing a time bin
bin_stc = nan(2*size(stc,1),4);
bin_stc(1,1:4) = stc(2,1:4);

current_tbin = 1; % refers to bin_stc row being filled
stc_last = 2;

dstc = diff(stc(:,1));
cstc = find(dstc ~= 0) + 1;
cstc(1:end-1,2) = cstc(2:end,1) - 1;
cstc(end,2) = size(stc,1);
cstc_track = 1;

while bin_stc(current_tbin, 1) + tbin_size <= stc(end,1)
    if ~conditions(stc_last) % do not allow any bin to include any ~condition stc rows
        bin_stc(current_tbin, 1) = stc(stc_last+1,1);
        stc_last = stc_last + 1;
        continue
    else
        while bin_stc(current_tbin, 1) < stc(stc_last+1,1)
            while ~any(cstc(cstc_track,1):cstc(cstc_track,2) == stc_last)
                cstc_track = cstc_track + 1;
            end
            match_idx = cstc(cstc_track,1):cstc(cstc_track,2); % to account for multiple simultaneously occupied bins
            match_b_idx = current_tbin:current_tbin-1+length(match_idx);
            if length(match_idx) > 1
                bin_stc(match_b_idx, 1) = bin_stc(current_tbin, 1);
                for i = 1:length(match_idx)
                    bin_stc(current_tbin+i-1, 2:4) = stc(match_idx(i), 2:4);
                end
                current_tbin = current_tbin + length(match_idx) - 1; % update index of latest filled tbin
            else
                bin_stc(current_tbin, 2:4) = stc(stc_last, 2:4);
            end
            
            bin_stc(current_tbin+1, 1) = bin_stc(current_tbin, 1) + tbin_size;
            current_tbin = current_tbin + 1;
        end
        while stc_last+1 <= size(stc,1) && stc(stc_last+1,1) <= bin_stc(current_tbin, 1)
            stc_last = stc_last + 1;
        end
    end
end

bin_stc = bin_stc(~isnan(bin_stc(:,1)),:);
bin_stc(:,5) = zeros(size(bin_stc,1),1); % 5th col contains number of spike

% place/view bin sieving
rows_remove = [];
for k = 1:size(bin_stc,1)
    if ~any(place_bins_sieved == bin_stc(k,2)) || ~any(view_bins_sieved == bin_stc(k,4))
        rows_remove = [rows_remove k];
    end
end
bin_stc(rows_remove,:) = [];

bin_dstc = diff(bin_stc(:,1));
bin_cstc = find(bin_dstc ~= 0) + 1;
bin_cstc(1:end-1,2) = bin_cstc(2:end,1) - 1;
bin_cstc(end,2) = size(bin_stc,1);

spike_latest = 1;
for k = 1:size(bin_cstc,1)
    for sp = spike_latest:length(spiketimes)
        if spiketimes(sp) >= bin_stc(bin_cstc(k,1),1) + tbin_size
            break
        elseif spiketimes(sp) >= bin_stc(bin_cstc(k,1),1) && spiketimes(sp) < bin_stc(bin_cstc(k,1),1) + tbin_size
            bin_stc(bin_cstc(k,1):bin_cstc(k,2), 5) = bin_stc(bin_cstc(k,1):bin_cstc(k,2), 5) + 1;
            spike_latest = sp + 1;
        end
    end
end

% generate duration maps for place and view
dstc = diff([stc(1:end,1); pv.data.rplmaxtime]);
stc = stc(conditions,:);
dstc = dstc(conditions);
place_dur = zeros(1600,1);
view_dur = zeros(5122,1);
for k = 1:size(stc,1)
    place_bin = stc(k,2);
    view_bin = stc(k,4);
    if place_bin >= 1 && place_bin <= 1600  % valid place bins
        place_dur(place_bin) = place_dur(place_bin) + dstc(k);
    end
    if view_bin >= 1 && view_bin <= 5122  % valid view bins
        view_dur(view_bin) = view_dur(view_bin) + dstc(k);
    end
end

% get place and view bins that have > 0 duration
place_good_bins = find(place_dur > 0);
view_good_bins = find(view_dur > 0);


vmpvData = struct;

vmpvData.ThresVel = ThresVel;
vmpvData.UseMinObs = UseMinObs;

vmpvData.bin_stc = bin_stc;
vmpvData.tbin_size = tbin_size;

vmpvData.place_good_bins = place_good_bins;
vmpvData.view_good_bins = view_good_bins;

% vmpvData.place_dur = place_dur;
% vmpvData.view_dur = view_dur;

save('vmpvData.mat','vmpvData','-v7.3');

end