function [ise_out] = ise1_3(actual_image, shuffled_images, dim1, dim2)
    %ISE using 1D histogram
    tic; %show how long ise caculation will take
    
    % parameters to discretize maps
    bin_resolution = 0.005;
    
    % binning each datapoint
    actual_disc = floor(actual_image/bin_resolution)+1;
    actual_disc(isnan(actual_disc)) = 0;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    shuffled_disc(isnan(shuffled_disc)) = 0;
    
    combined = [actual_disc; shuffled_disc]; %(1+shuffle) x dim1*dim2
    
    %Image spatial entropy calculation
    ISE=[];
    for  i=1:size(combined,1)
        if max(combined(i,:))==0 %an empty image don't have ISE
            disp('an image consist of NaN');
            ise_out = 0;
            ISE= [ISE ; ise_out];
            continue
        end
        %an image being discrtized by not small enough bin_reolution is
        %excluded from ISE calculation
        if max(combined(i,:))==1 
            ise_out = 0;
            ISE= [ISE ; ise_out];
            disp('floor(combined(i,:))=0');%every intensity value rounded to zero
            continue
        end
        % empty index array
        n=logical(zeros(1,size(combined(i,:),2)));
        %index
        image=zeros(1,size(combined(i,:),2));
        binLocations = [1:1:max(combined(i,:))];
        [counts] = histcounts(combined(i,:),binLocations);
        %convert count to probability
        counts=counts/sum(counts);
        s=size(counts,2);
        %assign propability of how often does a intensity occur in this image
        for ii=1:size(counts,2)
            if ii ~= s
                temp=(binLocations(ii)<=combined(i,:)) + (combined(i,:)<binLocations(ii+1));
                ind=temp==2;
            else
                temp=(binLocations(ii)<=combined(i,:)) + (combined(i,:)<=binLocations(ii+1));
                ind=temp==2;
            end
            image(ind)=counts(ii);
        end
        image = image/sum(image);
        image = reshape(image, dim1, dim2);

        %   H(X)
        centre = image; %X
        p_X= reshape(centre,[],1); %linearlize the map
        p_X(p_X==0) = []; %exclude the probability of having no intensity (NaN)
        
        [self_entropy] = entropy(p_X);
        
        ise_out = self_entropy;
        
        ISE= [ISE ; ise_out];
    end

    ise_out=ISE;
    disp(['time taken to calculate for ise1_1: ' num2str(toc)]);
end

function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);
    entropy= sum(entropy);%output sum of all entropy
    
end
