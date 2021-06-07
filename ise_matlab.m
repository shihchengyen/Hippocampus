function [ise_out] = ise_matlab(actual_image, shuffled_images, dim1, dim2)
    tic;
%     https://www.mathworks.com/help/images/ref/entropy.html
%     actual_image= [rand(1,8211)];
%     shuffled_images= [rand(10,8211)];
%     dim1=51;
%     dim2=161;
%     % parameters to discretize maps
%     bin_resolution = 0.05; % 0.005 
% 
%     % binning each datapoint
%     actual_disc = floor(actual_image/bin_resolution)+1;
%     actual_disc(isnan(actual_disc)) = 0;
%     shuffled_disc = floor(shuffled_images/bin_resolution)+1;
%     shuffled_disc(isnan(shuffled_disc)) = 0;
%     combined_disc = [actual_disc; shuffled_disc];
    combined_disc = [actual_image; shuffled_images];
    ISE=[];
    for  n=1:size(combined_disc,1)
        r=reshape(combined_disc(n,:),dim1,dim2);
        e=entropy(r);
        ISE= [ISE ; e];
    end
%     % reshaping to 2d structure, stacked by shuffles for 3d result
%     combined_disc = reshape(combined_disc, size(combined_disc,1), dim1, dim2);
%     
    ise_out=ISE;
    disp(['time taken to calculate for ise_matlab: ' num2str(toc)]);
end
   
% For find the value of entropy for a signal 1-dimensional, can use:
% p = hist(signal length);
% entropy = -sum(p.*log2(p));