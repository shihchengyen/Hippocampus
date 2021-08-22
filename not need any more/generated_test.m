ise = zeros(32, 14);

ii=1;
actual = reshape(b1,1,[]);
sh=[];
                dim1 = 40;
                dim2 = 40;
                address = 'C:\Users\teddy\Downloads\Data\New folder';
                smooth = [address,'\generated_'];
                %for ise1pixel
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);           
                    
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise1pixel_center.mat'];
                    end
                    save(filename,'ise'); 
                    c = c+2;
                end
ii=2;
actual = reshape(c1,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=3;
actual = reshape(c2,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=4;
actual = reshape(c3,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=5;
actual = reshape(c4,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=5;
actual = reshape(c4,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=6;
actual = reshape(c5,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=7;
actual = reshape(c6,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=8;
actual = reshape(c7,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=9;
actual = reshape(c8,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=10;
actual = reshape(c9,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=11;
actual = reshape(c10,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=12;
actual = reshape(c11,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=13;
actual = reshape(c12,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=14;
actual = reshape(c13,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=15;
actual = reshape(c14,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=16;
actual = reshape(c15,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=17;
actual = reshape(c16,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=18;
actual = reshape(c17,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=19;
actual = reshape(c18,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=20;
actual = reshape(c19,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=21;
actual = reshape(c20,1,[]);
sh=[];
dim1 = 33;
dim2 = 160;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=22;
actual = reshape(c21,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=23;
actual = reshape(c22,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=24;
actual = reshape(c23,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=25;
actual = reshape(c24,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=26;
actual = reshape(c25,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=27;
actual = reshape(c26,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth '.mat'];
end ise1pixel_center
save(filename,'ise');
c = c+2;
end

ii=28;
actual = reshape(c27,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=29;
actual = reshape(c28,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

 ii=30;
actual = reshape(c29,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end

ii=31;
actual = reshape(c30,1,[]);
sh=[];
dim1 = 40;
dim2 = 40;
address = 'C:\Users\teddy\Downloads\Data\New folder';
smooth = [address,'\generated_'];
%for ise1pixel
c=1;
bin=[0.1; 0.05; 0.01; 0.005];
for iii = 1: size(bin,1)
if iii>1
load(filename)
end
ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
ise(ii, c:c+1) = ise_out(1,:);
%save new .mat file
if iii==1
filename = [smooth 'ise1pixel_center.mat'];
end
save(filename,'ise');
c = c+2;
end
end
% % ii=32;
% % actual = reshape(b3,1,[]);
% % sh=[];
% %                 dim1 = 40;
% %                 dim2 = 40;
% %                 address = 'C:\Users\teddy\Downloads\Data\New folder';
% %                 smooth = [address,'\generated_'];
% %                 %for ise1pixel
% %                 c=1;
% %                 bin=[0.1; 0.05; 0.01; 0.005];
% %                 for iii = 1: size(bin,1)
% %                     if iii>1
% %                         load(filename)
% %                     end
% %                     ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
% %                     ise(ii, c:c+1) = ise_out(1,:);           
% %                     
% %                     %save new .mat file
% %                     if iii==1
% %                         filename = [smooth 'ise1pixel_center.mat'];
% %                     end
% %                     save(filename,'ise'); 
% %                     c = c+2;
% %                 end