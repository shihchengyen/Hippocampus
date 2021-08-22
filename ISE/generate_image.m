%Generated image
%Examine the result of vert_entropy, pos_angled_entropy,neg_angled_entropy,
%hor_entropy, self_entropy.
%result can be seen in 'test generated matrix's entropy' sheet
temp = rand(40,50); %random  number 0~1  %column c
% temp = randi([0 1],40,50); %random interger %column d
% temp=ones(40,50); %column e
% temp = [0 1 0 1;0 1 0 1;0 1 0 1;0 1 0 1;]; %column f
% temp=eye(4,5); %column g
% 
% % n = 4; %column h
% % e = ones(n,1);
% % A = spdiags([e -2*e e],-1:1,n,n);
% % temp=full(A);
% 
% % temp =[2 0 1 0 0; 0 2 0 1 0; 0 0 2 0 1;0 0 0 2 0]; %column i
% 
dim1=size(temp,1);
dim2=size(temp,2);
actual_image = reshape(temp,1,[]);
shuffled_images = [];
[vert_entropy, pos_angled_entropy,neg_angled_entropy, hor_entropy, self_entropy] = ise7(actual_image, shuffled_images, dim1, dim2);

% % %I did use the parts below to test how sensitive is the ise to the place
% % %field. But it turn out to be not that useful info. So can be ignore.
% % %place field--gradient 
% % %square
% % %%%%%%%%fix
% % place=ones(5,5)*0.1;
% % place=ones(4,4)*0.2;
% % place=ones(3,3)*0.3;
% % place=ones(2,2)*0.4;
% % %rectangular
% % place=ones(4,5)*0.1;
% % place=ones(3,4)*0.2;
% % place=ones(3,3)*0.3;
% % place=ones(2,2)*0.4;
% % %various diagnal
% % dia=[ 1 0 0;
% %     0 1 0;
% %     0 0 1];
% % dia=[ 3 2 1;
% %     2 3 2;
% %     1 2 3];
% % dia=[ 3 2 0;
% %     2 3 2;
% %     0 2 3];
% % dia(dia==0)=NaN;
% % 
% % %hole
% % hole = NaN(8,8);
% % %background
% % %sqaure 40 by 40
% % b1 = rand(40,40);
% % %rectangle 33 by 160
% % b2 = rand(33,160);
% % b3 = rand(40,40);
% % b3(8:15,8:15)= hole;
% % b3(24:31,8:15) = hole;
% % b3(8:15,24:31) = hole;
% % b3(24:31,24:31) = hole;
% % %place field--uniform square
% % pus1=ones(2,2)*0.9;
% % pus2=ones(3,3)*0.9;
% % pus3=ones(4,4)*0.9;
% % pus4=ones(5,5)*0.9;
% % pus5=ones(6,6)*0.9;
% % %place field--uniform rectangle
% % pus6=ones(2,3)*0.9;
% % pus7=ones(3,4)*0.9;
% % pus8=ones(4,5)*0.9;
% % pus9=ones(5,6)*0.9;
% % %place field--gradient square
% % pgs = ones(7,7)*0.5;
% % pgs(2:6,2:6) = ones(5,5)*0.6;
% % pgs(3:5,3:5) = ones(3,3)*0.9;
% % %center
% % c1 = b1; c1(20:21,20:21)=pus1;
% % c2 = b1; c2(20:22,20:22)=pus2;
% % c3 = b1; c3(20:23,20:23)=pus3;
% % c4 = b1; c4(20:24,20:24)=pus4;
% % c5 = b1; c5(20:25,20:25)=pus5;
% % c6 = b1; c6(20:21,20:22)=pus6;
% % c7 = b1; c7(20:22,20:23)=pus7;
% % c8 = b1; c8(20:23,20:24)=pus8;
% % c9 = b1; c9(20:24,20:25)=pus9;
% % c10 = b1; c10(20:26,20:26)=pgs;
% % 
% % c11 = b2; c11(20:21,20:21)=pus1;
% % c12 = b2; c12(20:22,20:22)=pus2;
% % c13 = b2; c13(20:23,20:23)=pus3;
% % c14 = b2; c14(20:24,20:24)=pus4;
% % c15 = b2; c15(20:25,20:25)=pus5;
% % c16 = b2; c16(20:21,20:22)=pus6;
% % c17 = b2; c17(20:22,20:23)=pus7;
% % c18 = b2; c18(20:23,20:24)=pus8;
% % c19 = b2; c19(20:24,20:25)=pus9;
% % c20 = b2; c20(20:26,20:26)=pgs;
% % 
% % c21 = b3; c21(20:21,20:21)=pus1;
% % c22 = b3; c22(20:22,20:22)=pus2;
% % c23 = b3; c23(20:23,20:23)=pus3;
% % c24 = b3; c24(20:24,20:24)=pus4;
% % c25 = b3; c25(20:25,20:25)=pus5;
% % c26 = b3; c26(20:21,20:22)=pus6;
% % c27 = b3; c27(20:22,20:23)=pus7;
% % c28 = b3; c28(20:23,20:24)=pus8;
% % c29 = b3; c29(20:24,20:25)=pus9;
% % c30 = b3; c30(20:26,20:26)=pgs;
% % %corner col
% % 
% % col1 = b1; col1(1:2,1:2)=pus1;
% % col2 = b1; col2(1:3,1:3)=pus2;
% % col3 = b1; col3(1:4,1:4)=pus3;
% % col4 = b1; col4(1:5,1:5)=pus4;
% % col5 = b1; col5(1:6,1:6)=pus5;
% % col6 = b1; col6(1:2,1:3)=pus6;
% % col7 = b1; col7(1:3,1:4)=pus7;
% % col8 = b1; col8(1:4,1:5)=pus8;
% % col9 = b1; col9(1:5,1:6)=pus9;
% % col10 = b1; col10(1:7,1:7)=pgs;
% % 
% % col11 = b2; col11(1:2,1:2)=pus1;
% % col12 = b2; col12(1:3,1:3)=pus2;
% % col13 = b2; col13(1:4,1:4)=pus3;
% % col14 = b2; col14(1:5,1:5)=pus4;
% % col15 = b2; col15(1:6,1:6)=pus5;
% % col16 = b2; col16(1:2,1:3)=pus6;
% % col17 = b2; col17(1:3,1:4)=pus7;
% % col18 = b2; col18(1:4,1:5)=pus8;
% % col19 = b2; col19(1:5,1:6)=pus9;
% % col20 = b2; col20(1:7,1:7)=pgs;
% % 
% % col21 = b3; col21(1:2,1:2)=pus1;
% % col22 = b3; col22(1:3,1:3)=pus2;
% % col23 = b3; col23(1:4,1:4)=pus3;
% % col24 = b3; col24(1:5,1:5)=pus4;
% % col25 = b3; col25(1:6,1:6)=pus5;
% % col26 = b3; col26(1:2,1:3)=pus6;
% % col27 = b3; col27(1:3,1:4)=pus7;
% % col28 = b3; col28(1:4,1:5)=pus8;
% % col29 = b3; col29(1:5,1:6)=pus9;
% % col30 = b3; col30(1:7,1:7)=pgs;
% % 
% % name = [b1; b2; b3; c1; c2; c3; c4; c5];
% % % %loop
% % % for n=1:size(name,1)
% % % %initialize the array for the month
% % %         
% % % %         shu = 10000; %shuffle images number
% % % %         ise_sh = NaN(shu,14,15); %(shuffle images number, 7bin*2, 15 cells)
% % % %         ise_2_5 = zeros(15, 14);
% % % %         ise_97 = zeros(15, 14);
% % % %         z =zeros(15, 14); %z-score
% % % %         h = zeros(15, 14); %lillietest
% % % %for smoothed
% % % %                 origin = vmp.data.origin;
% % %                 
% % % %                 sh=vmp.data.maps_adsmsh;
% % % %                 sh(sh==0)=NaN;
% % % ise = zeros(31, 14);
% % % ii=1;
% % % actual = reshape(b1,1,[]);
% % % sh=[];
% % %                 dim1 = 40;
% % %                 dim2 = 40;
% % %                 address = 'C:\Users\teddy\Downloads\Data\New folder';
% % %                 smooth = [address,'\generated_'];
% % %                 %for ise_QMRF
% % %                 c=1;
% % %                 bin=[0.1; 0.05; 0.01; 0.005];
% % %                 for iii = 1: size(bin,1)
% % %                     if iii>1
% % %                         load(filename)
% % %                     end
% % %                     ise_out = ise_QMRF(actual, sh, dim1, dim2, bin(iii));
% % %                     ise(ii, c:c+1) = ise_out(1,:);           
% % %                     
% % %                     %save new .mat file
% % %                     if iii==1
% % %                         filename = [smooth 'ise_QMRF_center.mat'];
% % %                     end
% % %                     save(filename,'ise'); 
% % %                     c = c+2;
% % %                 end
% % % end
% % % 
% % % b1 = rand(40,80);
% % % c1 = b1; c1(20:21,20:21)=pus1;
% % % c2 = b1; c2(20:22,20:22)=pus2;
% % % c3 = b1; c3(20:23,20:23)=pus3;
% % % c4 = b1; c4(20:24,20:24)=pus4;
% % % c5 = b1; c5(20:25,20:25)=pus5;
% % % c6 = b1; c6(20:21,20:22)=pus6;
% % % c7 = b1; c7(20:22,20:23)=pus7;
% % % c8 = b1; c8(20:23,20:24)=pus8;
% % % c9 = b1; c9(20:24,20:25)=pus9;
% % % c10 = b1; c10(20:26,20:26)=pgs;