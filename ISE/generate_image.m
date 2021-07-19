%fake image
% temp = rand(40,50); %random  number 0~1
% temp = randi([0 1],40,50); %random interger
% temp=ones(40,50);
% temp = [0 1 0 1;0 1 0 1;0 1 0 1;0 1 0 1;];
% temp=eye(4,5);

% n = 4;
% e = ones(n,1);
% A = spdiags([e -2*e e],-1:1,n,n);
% temp=full(A);

temp =[2 0 1 0 0; 0 2 0 1 0; 0 0 2 0 1;0 0 0 2 0];

dim1=size(temp,1);
dim2=size(temp,2);
actual_image = reshape(temp,1,[]);
shuffled_images = [];
% [vert_entropy, pos_angled_entropy,neg_angled_entropy, hor_entropy, self_entropy] = ise7(actual_image, shuffled_images, dim1, dim2);
[vert_entropy, pos_angled_entropy,neg_angled_entropy, hor_entropy, self_entropy] = ise8(actual_image, shuffled_images, dim1, dim2);
% [vert_entropy, pos_angled_entropy,neg_angled_entropy, hor_entropy, self_entropy] = ise9(actual_image, shuffled_images, dim1, dim2);
% [vert_entropy, pos_angled_entropy,neg_angled_entropy, hor_entropy, self_entropy] = ise10(actual_image, shuffled_images, dim1, dim2);



% dim1=size(temp,1);
% dim2=size(temp,2);
% actual_image = reshape(temp,1,[]);
% shuffled_images = [];
% [vert_entropy, pos_angled_entropy,neg_angled_entropy, hor_entropy, self_entropy] = ise7(actual_image, shuffled_images, dim1, dim2);
