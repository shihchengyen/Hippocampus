to_remove = '/hpctmp/e0407619';
to_replace = '/Volumes/Hippocampus';

cd /Volumes/Hippocampus/Data/picasso-misc
% all_vmpv = dir('**/vmpv.mat');
% disp('0')
% all_1vmpv = dir('**/1vmpv.mat');
% disp('1')
% all_100vmpv = dir('**/100vmpv.mat');
% disp('100')
% all_500vmpv = dir('**/500vmpv.mat');
% disp('500')

% for i = 1:length(all_vmpv)
%     cd(all_vmpv(i).folder);
%     load(all_vmpv(i).name);
%     old_string = pv.data.origin{1};
%     if contains(old_string, to_remove)
%         disp(old_string)
%         new_string = to_replace + extractBetween(old_string, length(to_remove) + 1, length(old_string));
%         disp(new_string)
%         pv.data.origin{1} = new_string;
%         save('vmpv.mat', 'pv', '-v7.3')
%     end
% end
% 
% for i = 1:length(all_1vmpv)
%     cd(all_1vmpv(i).folder);
%     load(all_1vmpv(i).name);
%     old_string = pv.data.origin{1};
%     if contains(old_string, to_remove)
%         disp(old_string)
%         new_string = to_replace + extractBetween(old_string, length(to_remove) + 1, length(old_string));
%         disp(new_string)
%         pv.data.origin{1} = new_string;
%         save('1vmpv.mat', 'pv', '-v7.3')
%     end
% end
% 
% for i = 1:length(all_100vmpv)
%     cd(all_100vmpv(i).folder);
%     load(all_100vmpv(i).name);
%     old_string = pv.data.origin{1};
%     if contains(old_string, to_remove)
%         disp(old_string)
%         new_string = to_replace + extractBetween(old_string, length(to_remove) + 1, length(old_string));
%         disp(new_string)
%         pv.data.origin{1} = new_string;
%         save('100vmpv.mat', 'pv', '-v7.3')
%     end
% end
% 
% for i = 1:length(all_500vmpv)
%     cd(all_500vmpv(i).folder);
%     load(all_500vmpv(i).name);
%     old_string = pv.data.origin{1};
%     if contains(old_string, to_remove)
%         disp(old_string)
%         new_string = to_replace + extractBetween(old_string, length(to_remove) + 1, length(old_string));
%         disp(new_string)
%         pv.data.origin{1} = new_string;
%         save('500vmpv.mat', 'pv', '-v7.3')
%     end
% end

cd /Volumes/Hippocampus/Data/picasso-misc
all_vmpc = dir('**/vmpc.mat');

for i = 1:length(all_vmpc)
    cd(all_vmpc(i).folder);
    load(all_vmpc(i).name);
    old_string = vmp.data.origin{1};
    if contains(old_string, to_remove)
        disp(old_string)
        new_string = to_replace + extractBetween(old_string, length(to_remove) + 1, length(old_string));
        disp(new_string)
        vmp.data.origin{1} = new_string;
        save('vmpc.mat', 'vmp', '-v7.3')
    end
end

cd /Volumes/Hippocampus/Data/picasso-misc
all_vmsv = dir('**/vmsv.mat');

for i = 1:length(all_vmsv)
    cd(all_vmsv(i).folder);
    load(all_vmsv(i).name);
    old_string = vms.data.origin{1};
    if contains(old_string, to_remove)
        disp(old_string)
        new_string = to_replace + extractBetween(old_string, length(to_remove) + 1, length(old_string));
        disp(new_string)
        vms.data.origin{1} = new_string;
        save('vmsv.mat', 'vms', '-v7.3')
    end
end