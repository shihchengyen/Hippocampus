% This function takes in the sorted spike timings 'firings.mda' and the
% curated cluster tags 'sorting-curation.json', applies the curation to the
% sorted spike timings, and returns the curated sorting results
% 'firings.curated.mda'. 
%
% This function is to be called from the channel
% directory of the sorting results (mountains/channel___), and the json
% file containing the curation labels is assumed to be in the Downloads
% folder (~/Downloads).

function curate_firings(exclude_labels)

    if ~exist('exclude_labels','var')
        exclude_labels = {'reject', 'noise', 'artifact'};
    end
    
    % Move 'sorting-curation.json' to the channel directory
    movefile ~/Downloads/sorting-curation.json .

    % Read in sorting and curation tags
    sorting = readmda('output/firings.mda');
    curation_tags = jsondecode(fileread('sorting-curation.json'));
    num_clusters = max(sorting(3,:));
    
    % Map new cluster labels to merged clusters
    new_cluster_labels = 1:num_clusters;
    for i = 1:size(curation_tags.mergeGroups,1)
        merge_group = curation_tags.mergeGroups(i,:);
        if iscell(merge_group)
            merge_group = cell2mat(merge_group);
        end
        for j = 1:length(merge_group)
            new_cluster_labels(merge_group(j)) = merge_group(1);
        end
    end
    
    % Look through all spike timings, tag spikes to be discarded and apply
    % merges
    for i = 1:size(sorting,2)
        clust = sorting(3,i);
        if ismember(curation_tags.labelsByUnit.(sprintf('x%d',clust)), exclude_labels)
           sorting(3,i) = -1;
        else
            sorting(3,i) = new_cluster_labels(clust);
        end
    end
    
    % Discard spikes
    discarded = (sorting(3,:) == -1);
    sorting(:, discarded) = [];
    
    % Export curated firings to .mda file
    writemda(sorting, 'firings.curated.mda', 'int32');
    
end
