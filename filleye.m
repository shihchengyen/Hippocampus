% This function takes in the eyelink messages (nx1) that contain text
% such as 'Cue Offset 24' or 'End Trial 35', along with eltimes (nx1) that
% contain the corresponding eyelink timestamps. The function fills in
% missing markers by comparing with ripple markers from the rpl object.
%
% Outputs elTrials, which is the filled in, reshaped timestamps (nx3).
% Outputs missing, which has the same (nx3) shape, filled with zeros except
% for markers that were originally not in eyelink.
% Outputs newMessages, which is a flattened, cell version of missing, with
% markers converted to text (eg. 34 now 'End Trial 34').
%
% The outputs for the function is as such, so as to make it compatible with
% the original implementation within completeData. This replacement handles
% all cases for missing eyelink markers, including multiple missing rows,
% missing markers from consecutive similar rows, which were cases where the
% previous implementation would have failed.
%
% used in completeData (~line 138)

function [elTrials, missing, newMessages] = filleye(messages, eltimes, rpl)


    eyelink_raw = NaN(1,length(messages));

    for i = 1:length(messages)
        full_text = messages{i};
        full_text = strsplit(full_text, ' ');
        full_text = full_text{length(full_text)}; 
        eyelink_raw(i) = str2double(full_text);
    end
    
    eye_timestamps = double(eltimes)';
    truth_timestamps = rpl.data.timeStamps;
    truth_timestamps = truth_timestamps * 1000;
    truth = rpl.data.markers;
    
    %%%%%% data transformed into format used by function %%%%%%%

    split_by_ones = NaN(2000,10);
    row = 1;
    col = 1;
    max_col = 1;
    start = 1;
    for i = 1:length(eyelink_raw) % naively splits the sequence by plausible triples or by cue onset
        if (eyelink_raw(i) < 20 && eyelink_raw(i) > 9) || col > 3
            row = row + 1;
            if col > max_col
                max_col = col;
            end
            col = 1;
        elseif (col~=1)
            if base < 20
                if eyelink_raw(i) ~= base + 10
                    if eyelink_raw(i) ~= base + 20
                        if eyelink_raw(i) ~= base + 30
                            row = row + 1;  
                            col = 1;
                        end
                    end
                end
            elseif base < 30
                if eyelink_raw(i) ~= base + 10
                    if eyelink_raw(i) ~= base + 20
                        row = row + 1;
                        col = 1;
                    end
                end
            else
                row = row + 1;
                col = 1;
            end    
        end
        split_by_ones(row, col) = eyelink_raw(i);
        base = eyelink_raw(i);
        col = col + 1;
        if start == 1
            start = 0;
        end
    end

    if sum(~isnan(split_by_ones(1,:))) ~= 0
        split_by_ones = split_by_ones(1:row,1:max_col-1);
    else
        split_by_ones = split_by_ones(2:row,1:max_col-1);
    end

    arranged_array = NaN(size(split_by_ones));

    for row = 1:size(split_by_ones,1)
        for col = 1:3
            if isnan(split_by_ones(row,col))
                break;
            end
            if split_by_ones(row,col) < 20
                arranged_array(row,1) = split_by_ones(row,col);
            elseif split_by_ones(row,col) < 30
                arranged_array(row,2) = split_by_ones(row,col);
            else
                arranged_array(row,3) = split_by_ones(row,col);
            end
        end
    end

    missing_rows = size(truth,1) - size(arranged_array,1);

    slice_after = NaN(missing_rows,2); % this section accounts for triples that look ok, but are made of two trials with the same posters
    slice_index = 1;
    for row = 1:size(arranged_array,1)
        if row > 316
            disp('debugger');
        end
        if ~isnan(arranged_array(row,1))
            if ~isnan(arranged_array(row,2))
                tmp = arranged_array';
                tmp = tmp(:);
                idx = sum(~isnan(tmp(1:3*(row-1)+1)));
                td = eye_timestamps(idx+1) - eye_timestamps(idx);

                    rpl_chunk = truth_timestamps(row:min([max([row row+missing_rows+1]) size(truth_timestamps,1)]),1:2);
                    rpl_chunk_flag = truth(row:min([max([row row+missing_rows+1]) size(truth_timestamps)]),1:2); %min was bug?

                rpl_chunk = rpl_chunk(rpl_chunk_flag(:,1)==arranged_array(row,1),:);           
                rpl_td = rpl_chunk(:,2) - rpl_chunk(:,1);

                if min(abs(rpl_td-td)) > 1500
                    slice_after(slice_index,:) = [row, 1];
                    slice_index = slice_index + 1;
                end
            elseif ~isnan(arranged_array(row,3))            
                tmp = arranged_array';
                tmp = tmp(:);
                idx = sum(~isnan(tmp(1:3*(row-1)+1))); 
                idx3 = sum(~isnan(tmp(1:3*(row-1)+3))); 
                td = eye_timestamps(idx3) - eye_timestamps(idx);

                    rpl_chunk = truth_timestamps(row:min([max([row row+missing_rows+1]) size(truth_timestamps,1)]),1:3);
                    rpl_chunk_flag = truth(row:min([max([row row+missing_rows+1]) size(truth_timestamps,1)]),1:3); %min was bug?

                rpl_chunk = rpl_chunk(rpl_chunk_flag(:,1)==arranged_array(row,1),:);  
                rpl_td = rpl_chunk(:,3) - rpl_chunk(:,1);

                if min(abs(rpl_td-td)) > 1500
                    slice_after(slice_index,:) = [row, 1];
                    slice_index = slice_index + 1;
                end            
            end
        elseif ~isnan(arranged_array(row,2))
            if ~isnan(arranged_array(row,3))
                tmp = arranged_array';
                tmp = tmp(:);
                idx = sum(~isnan(tmp(1:3*(row-1)+2))); 
                td = eye_timestamps(idx+1) - eye_timestamps(idx);

                    rpl_chunk = truth_timestamps(row:min([max([row row+missing_rows+1]) size(truth_timestamps,1)]),2:3);
                    rpl_chunk_flag = truth(row:min([max([row row+missing_rows+1]) size(truth_timestamps,1)]),2:3); %min was bug?

                rpl_chunk = rpl_chunk(rpl_chunk_flag(:,2)==arranged_array(row,2),:);
                rpl_td = rpl_chunk(:,2) - rpl_chunk(:,1);

                if min(abs(rpl_td-td)) > 1500
                    slice_after(slice_index,:) = [row, 1];
                    slice_index = slice_index + 1;
                end            
            end        
        end
    end

    slice_after = slice_after(1:slice_index-1,:);
    arranged_array = [arranged_array; NaN(missing_rows,3)];
    
%     disp('pre slicing');
%     arranged_array

    for slice = size(slice_after,1):-1:1 % slices according to previously identified segments
        new_array = NaN(size(arranged_array));
        new_array(1:slice_after(slice,1)-1,:) = arranged_array(1:slice_after(slice,1)-1,:);
        new_array(slice_after(slice,1),1:slice_after(slice,2)) = arranged_array(slice_after(slice,1),1:slice_after(slice,2));
        new_array(slice_after(slice,1)+1,slice_after(slice,2)+1:3) = arranged_array(slice_after(slice,1),slice_after(slice,2)+1:3);
        arranged_array(slice_after(slice,1)+1:size(arranged_array,1),:)
        new_array(slice_after(slice,1)+2:end,:) = arranged_array(slice_after(slice,1)+1:size(arranged_array,1)-1,:);
        arranged_array = new_array;
        missing_rows = missing_rows - 1;
    end

%     disp('post slicing');
%     arranged_array

    for row = 1:missing_rows % this segment attempts to identify where entire trials may have gone missing, by comparing with rpl timings
        error = nansum(truth - arranged_array,2);
        error_index = min(find(error~=0)); % insert before this
        if sum(abs(error)) == 0
            break;
        end
        if error_index == 1
            arranged_array = [NaN(1,3); arranged_array(1:end-1,:)];
        else
            for col = 1:3
                if ~isnan(arranged_array(error_index-1,col))
                    pre_id = rem(arranged_array(error_index-1,col),10);
                    break;
                end
            end % identify of the preceeding trial determined
            % looking up how many trials before this have the same identity
            count = 0;
            while 1
                if error_index-1-count == 0
                    break;
                end
                for col2 = 1:3
                    if ~isnan(arranged_array(error_index-1-count,col2))
                        pre_id_check = rem(arranged_array(error_index-1-count,col2),10);
                        break;
                    end
                end
                if pre_id_check ~= pre_id
                    break;
                end
                if error_index-2-count == 0
                    break;
                end
                count = count + 1;
            end
            % count now stores the number of repeated posters before the
            % misalignment has been detected (need to test all possible
            % locations).
            disp(count);

            eye_start_trials = NaN(count+2,1);
            eye_start_count = 1;

            esi = 0;
            for r = 1:size(arranged_array,1)
                for c = 1:3
                    if ~isnan(arranged_array(r,c))
                        esi = esi + 1;
                    end
                    if (r >= error_index-count-1) && (r <= error_index)
                        if c == 1
                            if ~isnan(arranged_array(r,c))
                                eye_start_trials(eye_start_count,1) = eye_timestamps(esi);
                            elseif ~isnan(arranged_array(r,c+1))
                                disp('taking cue offset and cutting 2seconds to estimate start trial timing');
                                eye_start_trials(eye_start_count,1) = eye_timestamps(esi+1)-2000;
                            else
                                disp('taking end trial and cutting 10seconds to estimate start trial timing');
                                eye_start_trials(eye_start_count,1) = eye_timestamps(esi+1)-10000;
                            end
                            eye_start_count = eye_start_count + 1;
                        end
                    end
                end
            end

                rpl_start_trials = truth_timestamps(error_index-count-1:error_index,1);
                diff_eye = diff(eye_start_trials);
                diff_rpl = diff(rpl_start_trials);
                discrepency = diff_eye - diff_rpl;
                [~,row_to_insert] = max(discrepency);
                
                arranged_array = [arranged_array(1:error_index-count-2+row_to_insert,:); NaN(1,3); arranged_array(error_index-count-1+row_to_insert:end,:)];
                arranged_array = arranged_array(1:end-1,:);
        end
    end

    if nansum(abs(double(arranged_array) - double(truth))) > 0
        clear error;
        error('eyelink was not properly arranged. current arrangement still clashes with ripple')
    end
    
    missing = truth.*double(isnan(arranged_array)); %%%%% ready for output
    newMessages = cell(3*size(truth,1),1);
    flat_truth = truth';
    flat_truth = flat_truth(:);
    flat_truth_time = truth_timestamps';
    flat_truth_time = flat_truth_time(:);
    flat_eye = arranged_array';
    flat_eye = flat_eye(:);
    flat_truth = flat_truth.*double(isnan(flat_eye));
    for i = 1:length(flat_truth)
        if flat_truth(i) ~= 0
            if flat_truth(i) < 20
                text = ['Start Trial ' num2str(flat_truth(i))];
                newMessages{i} = text;
            elseif flat_truth(i) < 30
                text = ['Cue Offset ' num2str(flat_truth(i))];
                newMessages{i} = text;                
            elseif flat_truth(i) < 40
                text = ['End Trial ' num2str(flat_truth(i))];
                newMessages{i} = text;                
            else
                text = ['Timeout ' num2str(flat_truth(i))];
                newMessages{i} = text;                
            end
        end
    end %%%%% ready for output

    elTrials = NaN(1,3*size(missing,1));
    counter = 1;
    for i = 1:length(flat_eye)
        if ~isnan(flat_eye(i))
            elTrials(i) = eltimes(counter);
            counter = counter + 1;
        end
    end
    for i = 1:length(elTrials)
        if isnan(elTrials(i))
            if i == 1
                inv_delta = flat_truth_time(i+1) - flat_truth_time(i);
                elTrials(i) = round(elTrials(i+1) - inv_delta);
                disp('shouldnt see nans here');
                elTrials(i)
            else
                delta = flat_truth_time(i) - flat_truth_time(i-1);
                elTrials(i) = round(elTrials(i-1) + delta);
                disp('shouldnt see nans here');
                elTrials(i)
            end
        end
    end    
    elTrials = reshape(elTrials, 3, [])'; %%%%% ready for output

end






