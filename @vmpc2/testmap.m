function maps_adsm = testmap()

Args.GridSteps = 40;
Args.Alpha = 10000;
Args.NumShuffles = 0;

maps_raw_data = load('/Users/celinewang/Desktop/maps_raw.mat');
spikes_count_data = load('/Users/celinewang/Desktop/spikes_count.mat');
dur_raw_data = load('/Users/celinewang/Desktop/dur_raw.mat');

maps_raw = struct2array(maps_raw_data);
spikes_count = struct2array(spikes_count_data);
dur_raw = struct2array(dur_raw_data);

[maps_adsm, ~, ~, ~, ~, ~, ~, ~] = smoothMaps(maps_raw, dur_raw, spikes_count, spikes_count, Args);

figure;
imagesc(reshape(maps_adsm(1,:), [Args.GridSteps, Args.GridSteps]));
axis image;
title('maps\_adsm');
colorbar;

end
