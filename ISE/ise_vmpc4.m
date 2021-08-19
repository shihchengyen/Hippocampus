clear
clc

cwd ='C:\Users\Teddy\Downloads\data\New folder'; %ASUS
% cwd = 'C:\Users\teddy\Downloads\Data\New folder'; %hp
cd(cwd);
list=['20181031'; '20181101';'20181102'];
%For place map
dim1 = 40;
dim2 = 40;
for i=2:size(list,1)
    list=['20181031'; '20181101';'20181102'];
    
    cd(list(i,:));
    if i==1 %1031
        %restart smoothed map at ch30c1 next
%        list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
%        list1031=['ch35c1'];
       list1031=['ch19c1'];
       %initialize the array for the month
        ise = zeros(size(list1031,1), 120);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,120,size(list1031,1)); %(shuffle images number, version, 15 cells?)
        ise_2_5 = zeros(size(list1031,1), 120);
        ise_97 = zeros(size(list1031,1), 120);
        z =zeros(size(list1031,1), 120); %z-score
        h = zeros(size(list1031,1), 120); %lillietest

       for ii=1:size(list1031,1)
            cd(list1031(ii,:));
%             cwd ='C:\Users\Teddy\Downloads\data\New folder';
            cwd = 'C:\Users\teddy\Downloads\Data\New folder';
            address = [cwd, '\', list(i,:)];
            
            %ISE
                load('vmpc.mat');   
                
                %for smoothed
                origin = vmp.data.origin;
                actual = vmp.data.maps_adsm;
                l=isnan(actual);
                v=~l;
                sh=NaN(shu,dim1*dim2);
                sh(:,v)=vmp.data.maps_adsmsh(:,v);
                smooth = [address,'\smoothed_'];
                
                %for ise_threshold
                c=1;
                bin = [0.1; 0.05; 0.01; 0.005; 0];
                divide = [0; 10; 15; 20];
                threshold = [0; 90; 70; 50];
                change = [2; 1; 0];
                pixel = [1; 2];
                filename = [smooth 'ise_threshold_p_19c1_n.mat'];
                if ii>1
                        load(filename)
                end
                for iii = 1: size(bin,1)
                    for i4 = 1: size(divide,1)
                        for i5 = 1: size(threshold,1)
                            for i6 = 1: size(change,1)
                                for i7 = 1: size(pixel,1)
                                    if bin(iii)==0 && divide(i4)==0
                                        ise(ii, c) = NaN;
                                        ise_sh(:,c,ii) = NaN(shu,1);
                                        ise_2_5(ii, c) = NaN;
                                        ise_97(ii, c) = NaN;
                                        z(ii, c) = NaN;
                                        h(ii, c)=2;
                                        continue
                                    elseif bin(iii)~=0 && divide(i4)~=0
                                        ise(ii, c) = NaN;
                                        ise_sh(:,c,ii) = NaN(shu,1);
                                        ise_2_5(ii, c) = NaN;
                                        ise_97(ii, c) = NaN;
                                        z(ii, c) = NaN;
                                        h(ii, c)=2;
                                        continue   
                                    end
                                    ise_out = ise_threshold(actual, sh, dim1, dim2, bin(iii), divide(i4), threshold(i5), change(i6), pixel(i7));
                                    ise(ii, c) = ise_out(1);
                                    ise_sh(:,c,ii) = ise_out(2:end);

                                    ise_sh1 = ise_sh(:,c,ii);
                                    ise_sh1(ise_sh1==0) = []; %exclude zero ise

                                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                                    if sum(ise_sh1)>0
                                        h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                                    else
                                        h(ii, c)=2; %error, an empty array 
                                    end
                                    c = c+1;
                                end
                            end
                        end
                    end
                end
                 %save new .mat file
                 save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
            cd ..
       end
       cd .. %back one directory
    end
    
    
    if i==2 %1101
%        list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch29c4';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1'];
%        list1101=['ch19c1';'ch19c2';'ch21c1';'ch29c1';'ch29c2';'ch30c2'];
       list1101=['ch30c2'];
       %initialize the array for the month
        ise = zeros(size(list1101,1), 120);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,120,size(list1101,1)); %(shuffle images number, version, 15 cells?)
        ise_2_5 = zeros(size(list1101,1), 120);
        ise_97 = zeros(size(list1101,1), 120);
        z =zeros(size(list1101,1), 120); %z-score
        h = zeros(size(list1101,1), 120); %lillietest
        
       for ii=1:size(list1101,1)
           cd(list1101(ii,:));
           cwd ='C:\Users\Teddy\Downloads\data\New folder';
%            cwd = 'C:\Users\teddy\Downloads\Data\New folder';
           address = [cwd, '\', list(i,:)];
            
           %ISE
                load('vmpc.mat');   
                
                %for smoothed
                origin = vmp.data.origin;
                actual = vmp.data.maps_adsm;
                l=isnan(actual);
                v=~l;
                sh=NaN(shu,dim1*dim2);
                sh(:,v)=vmp.data.maps_adsmsh(:,v);
                smooth = [address,'\smoothed_'];
                
                %for ise_threshold
                c=1;
                bin = [0.1; 0.05; 0.01; 0.005; 0];
                divide = [0; 10; 15; 20];
                threshold = [0; 90; 70; 50];
                change = [2; 1; 0];
                pixel = [1; 2];
                filename = [smooth 'ise_threshold_p_30c2.mat'];
                if ii>1
                        load(filename)
                end
                for iii = 1: size(bin,1)
                    for i4 = 1: size(divide,1)
                        for i5 = 1: size(threshold,1)
                            for i6 = 1: size(change,1)
                                for i7 = 1: size(pixel,1)
                                    if bin(iii)==0 && divide(i4)==0
                                        ise(ii, c) = NaN;
                                        ise_sh(:,c,ii) = NaN(shu,1);
                                        ise_2_5(ii, c) = NaN;
                                        ise_97(ii, c) = NaN;
                                        z(ii, c) = NaN;
                                        h(ii, c)=2;
                                        continue
                                    elseif bin(iii)~=0 && divide(i4)~=0
                                        ise(ii, c) = NaN;
                                        ise_sh(:,c,ii) = NaN(shu,1);
                                        ise_2_5(ii, c) = NaN;
                                        ise_97(ii, c) = NaN;
                                        z(ii, c) = NaN;
                                        h(ii, c)=2;
                                        continue   
                                    end
                                    ise_out = ise_threshold(actual, sh, dim1, dim2, bin(iii), divide(i4), threshold(i5), change(i6), pixel(i7));
                                    ise(ii, c) = ise_out(1);
                                    ise_sh(:,c,ii) = ise_out(2:end);

                                    ise_sh1 = ise_sh(:,c,ii);
                                    ise_sh1(ise_sh1==0) = []; %exclude zero ise

                                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                                    if sum(ise_sh1)>0
                                        h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                                    else
                                        h(ii, c)=2; %error, an empty array 
                                    end

                                    c = c+1;
                                end
                            end
                        end
                    end
                end
                 %save new .mat file
                 save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
            cd ..
       end
       cd .. %back one directory
    end
    
    if i==3 %1102
%        list1102=['ch09c1';'ch19c1';'ch19c2';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch31c1';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
       list1102=['ch31c1'];
       %initialize the array for the month
        ise = zeros(size(list1102,1), 120);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,120,size(list1102,1)); %(shuffle images number, version, 15 cells?)
        ise_2_5 = zeros(size(list1102,1), 120);
        ise_97 = zeros(size(list1102,1), 120);
        z =zeros(size(list1102,1), 120); %z-score
        h = zeros(size(list1102,1), 120); %lillietest
        
       for ii=1:size(list1102,1)
           cd(list1102(ii,:));
          cwd ='C:\Users\Teddy\Downloads\data\New folder';
%            cwd = 'C:\Users\teddy\Downloads\Data\New folder';
            address = [cwd, '\', list(i,:)];
            
            %ISE
                load('vmpc.mat');   
                
                %for smoothed
                origin = vmp.data.origin;
                actual = vmp.data.maps_adsm;
                l=isnan(actual);
                v=~l;
                sh=NaN(shu,dim1*dim2);
                sh(:,v)=vmp.data.maps_adsmsh(:,v);
                smooth = [address,'\smoothed_'];
                
                %for ise_threshold
                c=1;
                bin = [0.1; 0.05; 0.01; 0.005; 0];
                divide = [0; 10; 15; 20];
                threshold = [0; 90; 70; 50];
                change = [2; 1; 0];
                pixel = [1; 2];
                filename = [smooth 'ise_threshold_p_31c1.mat'];
                if ii>1
                        load(filename)
                end
                for iii = 1: size(bin,1)
                    for i4 = 1: size(divide,1)
                        for i5 = 1: size(threshold,1)
                            for i6 = 1: size(change,1)
                                for i7 = 1: size(pixel,1)
                                    if bin(iii)==0 && divide(i4)==0
                                        ise(ii, c) = NaN;
                                        ise_sh(:,c,ii) = NaN(shu,1);
                                        ise_2_5(ii, c) = NaN;
                                        ise_97(ii, c) = NaN;
                                        z(ii, c) = NaN;
                                        h(ii, c)=2;
                                        continue
                                    elseif bin(iii)~=0 && divide(i4)~=0
                                        ise(ii, c) = NaN;
                                        ise_sh(:,c,ii) = NaN(shu,1);
                                        ise_2_5(ii, c) = NaN;
                                        ise_97(ii, c) = NaN;
                                        z(ii, c) = NaN;
                                        h(ii, c)=2;
                                        continue   
                                    end
                                    ise_out = ise_threshold(actual, sh, dim1, dim2, bin(iii), divide(i4), threshold(i5), change(i6), pixel(i7));
                                    ise(ii, c) = ise_out(1);
                                    ise_sh(:,c,ii) = ise_out(2:end);

                                    ise_sh1 = ise_sh(:,c,ii);
                                    ise_sh1(ise_sh1==0) = []; %exclude zero ise

                                    ise_2_5(ii, c) = prctile(ise_sh1, 2.5);
                                    ise_97(ii, c) = prctile(ise_sh1, 97.5);
                                    z(ii, c) =(ise(ii, c)-mean(ise_sh1))/std(ise_sh1);
                                    if sum(ise_sh1)>0
                                        h(ii, c)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                                    else
                                        h(ii, c)=2; %error, an empty array 
                                    end
                                    c = c+1;
                                end
                            end
                        end
                    end
                end
                 %save new .mat file
                 save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h'); 
            cd ..
       end
       cd .. %back one directory
    end
end
% % % map
% load('vmpc.mat')
% figure; plotmap(vmp.data.maps_raw,'place')
% figure; plotmap(vmp.data.maps_adsm,'place')
% figure; plotmap(vmp.data.maps_adsm1,'place')
% figure; plotmap(vmp.data.maps_adsm2,'place')