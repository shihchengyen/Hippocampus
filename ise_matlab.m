function [ise_out] = ise_matlab(actual_image, shuffled_images, dim1, dim2)
    tic;
%     https://www.mathworks.com/help/images/ref/entropy.html
%     for testing:
%     actual_image= [rand(1,8211)];
%     shuffled_images= [rand(10,8211)];
%     dim1=51;
%     dim2=161;
    % parameters to discretize maps
    bin_resolution = 0.05; % 0.005 

    % binning each datapoint
    actual_disc = floor(actual_image/bin_resolution)+1;
    actual_disc(isnan(actual_disc)) = 0;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    shuffled_disc(isnan(shuffled_disc)) = 0;
    
combined = [actual_disc; shuffled_disc];
[counts,binLocations] = imhist(actual_disc);
%  assign propability
%convert count to probability
counts=counts/sum(counts);
% empty index array
n=logical(zeros(size(combined,1),size(combined,2)));
%index
s=size(counts,1);
image=zeros(size(combined,1),size(combined,2));
for c=1:size(combined,1)
    i=1;
    for i=1:size(counts)
        if i ~= s
            temp=(binLocations(i)<=combined(c,:)) + (combined(c,:)<binLocations(i+1));
            ind=temp==2;
            n(c,:)=ind;
        else
            temp=(binLocations(i-1)<=combined(c,:)) + (combined(c,:)<=binLocations(i));
            ind=temp==2;
            n(c,:)=ind;
        end
        image(n(c,:))=counts(i);
    end
end

    ISE=[];
    for  i=1:size(combined,1)
        r=image(i,:);
        e=entropy(r);
        ISE= [ISE ; e];
    end

    ise_out=ISE;
    disp(['time taken to calculate for ise_matlab: ' num2str(toc)]);
end
   