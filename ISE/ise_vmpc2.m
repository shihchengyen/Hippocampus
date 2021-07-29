clear
clc

% cwd ='C:\Users\Teddy\Downloads\data\New folder'; %ASUS
cwd = 'C:\Users\teddy\Downloads\Data\New folder'; %hp
cd(cwd);
list=['20181031'; '20181101';'20181102'];
% % % % % % % % % % STOP AT 1102 ch09c1

for i=2:size(list,1)
    list=['20181031'; '20181101';'20181102'];
    
    cd(list(i,:));
    if i==1 %1031
        %restart smoothed map at ch30c1 next
       list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
       
       %initialize the array for the month
        ise = zeros(15, 14);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,14,15); %(shuffle images number, 7bin*2, 15 cells)
        ise_2_5 = zeros(15, 14);
        ise_97 = zeros(15, 14);
        z =zeros(15, 14); %z-score
        h = zeros(15, 14); %lillietest

       for ii=1:size(list1031,1)
            list=['20181031'; '20181101';'20181102'];
            list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
            cd(list1031(ii,:));
%             cwd ='C:\Users\Teddy\Downloads\data\New folder';
            cwd = 'C:\Users\teddy\Downloads\Data\New folder';
            address = [cwd, '\', list(i,:)];
            
            %ISE
                load('vmpc.mat');   
                
                %for smoothed
                origin = vmp.data.origin;
                actual = vmp.data.maps_adsm;
                sh=vmp.data.maps_adsmsh;
                sh(sh==0)=NaN;
                dim1 = 40;
                dim2 = 40;
                smooth = [address,'\smoothed_'];
                %for ise_QMRF
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise_QMRF(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise_QMRF_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise1pixel
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise1pixel_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise2pixels
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise2pixels(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise2pixels_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise3pixels
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise3pixels(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise3pixels_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
%                 %for raw
%%%%%%%%%%%%%%need to add
%                 load('vmpc_ted.mat')
%                 origin = vmp.data.origin;
%                 actual = vmp.data.maps_raw;
%                 sh=vmp.data.maps_rawsh;%not yet created
%                 sh(sh==0)=NaN;
%                 dim1 = 40;
%                 dim2 = 40;
%                 raw = [address,'\raw_'];
%                 %for ise_QMRF
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise_QMRF(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise_QMRF_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise1pixel
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise1pixel_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise2pixels
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise2pixels(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise2pixels_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise3pixels
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise3pixels(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise3pixels_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
            cd ..
       end
       cd .. %back one directory
    end
    
    if i==2 %1101
       list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch29c4';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1'];
       
       %initialize the array for the month
        ise = zeros(13, 14);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,14,13); %(shuffle images number, 7bin*2, 15 cells)
        ise_2_5 = zeros(13, 14);
        ise_97 = zeros(13, 14);
        z =zeros(13, 14); %z-score
        h = zeros(13, 14); %lillietest
       for ii=10:size(list1101,1)
           list=['20181031'; '20181101';'20181102'];
           list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch29c4';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1']; 
           cd(list1101(ii,:));
%            cwd ='C:\Users\Teddy\Downloads\data\New folder';
           cwd = 'C:\Users\teddy\Downloads\Data\New folder';
           address = [cwd, '\', list(i,:)];
            
           
            
            %ISE
                load('vmpc.mat');   
                
                %for smoothed
                origin = vmp.data.origin;
                actual = vmp.data.maps_adsm;
                sh=vmp.data.maps_adsmsh;
                sh(sh==0)=NaN;
                dim1 = 40;
                dim2 = 40;
                smooth = [address,'\smoothed_'];
                %for ise_QMRF
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise_QMRF(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise_QMRF_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise1pixel
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise1pixel_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise2pixels
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise2pixels(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise2pixels_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise3pixels
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise3pixels(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise3pixels_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
%                 %for raw
%                 load('vmpc_ted.mat')
%                 origin = vmp.data.origin;
%                 actual = vmp.data.maps_raw;
%                 sh=vmp.data.maps_rawsh;%not yet created
%                 sh(sh==0)=NaN;
%                 dim1 = 40;
%                 dim2 = 40;
%                 raw = [address,'\raw_'];
%                 %for ise_QMRF
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise_QMRF(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise_QMRF_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise1pixel
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise1pixel_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise2pixels
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise2pixels(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise2pixels_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise3pixels
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise3pixels(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise3pixels_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
            cd ..
       end
       cd .. %back one directory
    end
    
    if i==3 %1102
       list1102=['ch09c1';'ch19c1';'ch19c2';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch31c1';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
       %initialize the array for the month
        ise = zeros(13, 14);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,14,13); %(shuffle images number, 7bin*2, 15 cells)
        ise_2_5 = zeros(13, 14);
        ise_97 = zeros(13, 14);
        z =zeros(13, 14); %z-score
        h = zeros(13, 14); %lillietest
       for ii=1:size(list1102,1)
           list=['20181031'; '20181101';'20181102'];
           list1102=['ch09c1';'ch19c1';'ch19c2';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch31c1';'ch43c1';'ch43c2';'ch45c1';'ch45c2']; 
           cd(list1102(ii,:));
%           cwd ='C:\Users\Teddy\Downloads\data\New folder';
           cwd = 'C:\Users\teddy\Downloads\Data\New folder';
            address = [cwd, '\', list(i,:)];
            
           
            %ISE
                load('vmpc.mat');   
                
                %for smoothed
                origin = vmp.data.origin;
                actual = vmp.data.maps_adsm;
                sh=vmp.data.maps_adsmsh;
                sh(sh==0)=NaN;
                dim1 = 40;
                dim2 = 40;
                smooth = [address,'\smoothed_'];
                %for ise_QMRF
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise_QMRF(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise_QMRF_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise1pixel
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise1pixel_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise2pixels
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise2pixels(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise2pixels_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
                %for ise3pixels
                c=1;
                bin=[0.1; 0.05; 0.01; 0.005];
                for iii = 1: size(bin,1)
                    if iii>1
                        load(filename)
                    end
                    ise_out = ise3pixels(actual, sh, dim1, dim2, bin(iii));
                    ise(ii, c:c+1) = ise_out(1,:);
                    ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,c,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    ise_sh2 = ise_sh(:,c+1,ii);
                    ise_sh2(ise_sh2==0) = []; %exclude zero ise

                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                    ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                    ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                    z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
                    h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                    h(ii, c+1)=lillietest(ise_sh2);
                    %save new .mat file
                    if iii==1
                        filename = [smooth 'ise3pixels_p.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
                    c = c+2;
                end
                
%                 %for raw
%                 load('vmpc_ted.mat')
%                 origin = vmp.data.origin;
%                 actual = vmp.data.maps_raw;
%                 sh=vmp.data.maps_rawsh;%not yet created
%                 sh(sh==0)=NaN;
%                 dim1 = 40;
%                 dim2 = 40;
%                 raw = [address,'\raw_'];
%                 %for ise_QMRF
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise_QMRF(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise_QMRF_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise1pixel
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise1pixel(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise1pixel_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise2pixels
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise2pixels(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise2pixels_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
%                 
%                 %for ise3pixels
%                 c=1;
%                 bin=[0.1; 0.05; 0.01; 0.005];
%                 for iii = 1: size(bin,1)
%                     if iii>1
%                         load(filename)
%                     end
%                     ise_out = ise3pixels(actual, sh, dim1, dim2, bin(iii));
%                     ise(ii, c:c+1) = ise_out(1,:);
%                     ise_sh(:,c:c+1,ii) = ise_out(2:end,:);
%                    
%                     ise_sh1 = ise_sh(:,c,ii);
%                     ise_sh1(ise_sh1==0) = []; %exclude zero ise
%                     ise_sh2 = ise_sh(:,c+1,ii);
%                     ise_sh2(ise_sh2==0) = []; %exclude zero ise
% 
%                     ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
%                     ise_2_5(ii, c+1) = prctile(ise_sh2, 2.5);
%                     ise_97(ii, c) = prctile(ise_sh1, 97.5);
%                     ise_97(ii, c+1) = prctile(ise_sh2, 97.5);
%                     z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
%                     z(ii, c+1) =(ise(ii, c+1)-mean(ise_sh2))/std(ise_sh2);
%                     h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
%                     h(ii, c+1)=lillietest(ise_sh2);
%                     %save new .mat file
%                     if iii==1
%                         filename = [raw 'ise3pixels_p.mat'];
%                     end
%                     save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
%                     c = c+2;
%                 end
            cd ..
       end
       cd .. %back one directory
    end
end
