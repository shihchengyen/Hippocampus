days = [20180907
20180906
20180905];

error_days = [];

for i = 1:length(days)
    ori_path = pwd;
    cd(num2str(days(i)));
    try
        mountain_batch('hippo');
        cd('..');
    catch
        error_days = [error_days; days(i)];
        cd(ori_path);
    end
end

for i = 1:length(days)
    cd(num2str(days(i)));
    cd('mountains');
    [a, b] = unix('find . -name "channel???" | wc -l');
    [c, d] = unix('find . -name "firings.mda" | wc -l');
    [e, f] = unix('find . -name "metrics.json" | wc -l');
    disp(strcat(num2str(days(i)), b, d, f));
    cd('../..');
end

disp(['days with issues: ' num2str(length(error_days))]);
disp(error_days);


