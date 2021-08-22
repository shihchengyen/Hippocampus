cwd ='C:\Users\Teddy\Downloads\data\New folder'; %ASUS
% cwd = 'C:\Users\teddy\Downloads\Data\New folder'; %hp
cd(cwd);
list=['20181031'; '20181101';'20181102'];
dim1 = 40;
                dim2 = 40;
% % % % % % % % % % STOP AT 1102 ch09c1

for i=1:size(list,1)
    list=['20181031'; '20181101';'20181102'];
    
    cd(list(i,:));
    if i==1 %1031
       list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
       
       %initialize the array for the month
        ise = zeros(15, 1);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,1,15); %(shuffle images number, 7bin*2, 15 cells)
        ise_2_5 = zeros(15, 1);
        ise_97 = zeros(15, 1);
        z =zeros(15, 1); %z-score
        h = zeros(15, 1); %lillietest

       for ii=1:size(list1031,1)
            list=['20181031'; '20181101';'20181102'];
            list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
            cd(list1031(ii,:));
            cwd ='C:\Users\Teddy\Downloads\data\New folder';
%             cwd = 'C:\Users\teddy\Downloads\Data\New folder';
            address = [cwd, '\', list(i,:)];
            
            %ISE
                load('vmpc.mat');   
                
                %for smoothed
                origin = vmp.data.origin;
                actual = vmp.data.maps_adsm;
                l=isnan(actual);
                v=~l;
                sh=NaN(10000,dim1*dim2);
                sh(:,v)=vmp.data.maps_adsmsh(:,v);
                
                dim1 = 40;
                dim2 = 40;
                smooth = [address,'\smoothed_'];
                %for ise1pixel3
                    if ii>1
                        load(filename)
                    end
                    ise_out = ise1pixel3(actual, sh, dim1, dim2);
                    ise(ii) = ise_out(1,:);
                    ise_sh(:,1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,1,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    
                    ise_2_5(ii, 1) = prctile(ise_sh1, 2.5);
                    ise_97(ii, 1) = prctile(ise_sh1, 97.5);
                    
                    z(ii, 1) =(ise(ii)-mean(ise_sh1))/std(ise_sh1);
                    
                    h(ii, 1)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                   
                    %save new .mat file
                    if ii==1
                        filename = [smooth 'ise1pixel3_p50_15.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h');               
                
               cd ..
       end
       cd .. %back one directory
    end
    
    if i==2 %1101
       list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch29c4';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1'];
       
       %initialize the array for the month
        ise = zeros(13, 1);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,1,13); %(shuffle images number, 7bin*2, 15 cells)
        ise_2_5 = zeros(13, 1);
        ise_97 = zeros(13, 1);
        z =zeros(13, 1); %z-score
        h = zeros(13, 1); %lillietest
       for ii=1:size(list1101,1)
           list=['20181031'; '20181101';'20181102'];
           list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch29c4';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1']; 
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
                sh=NaN(10000,dim1*dim2);
                sh(:,v)=vmp.data.maps_adsmsh(:,v);
                
                dim1 = 40;
                dim2 = 40;
                smooth = [address,'\smoothed_'];
                %for ise1pixel3
                    if ii>1
                        load(filename)
                    end
                    ise_out = ise1pixel3(actual, sh, dim1, dim2);
                    ise(ii) = ise_out(1,:);
                    ise_sh(:,1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,1,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    
                    ise_2_5(ii, 1) = prctile(ise_sh1, 2.5);
                    ise_97(ii, 1) = prctile(ise_sh1, 97.5);
                    
                    z(ii, 1) =(ise(ii)-mean(ise_sh1))/std(ise_sh1);
                    
                    h(ii, 1)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                   
                    %save new .mat file
                    if ii==1
                        filename = [smooth 'ise1pixel3_p50_15.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h');
            cd ..
       end
       cd .. %back one directory
    end
    
    if i==3 %1102
       list1102=['ch09c1';'ch19c1';'ch19c2';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch31c1';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
       %initialize the array for the month
        ise = zeros(13, 1);
        shu = 10000; %shuffle images number
        ise_sh = NaN(shu,1,13); %(shuffle images number, 7bin*2, 15 cells)
        ise_2_5 = zeros(13, 1);
        ise_97 = zeros(13, 1);
        z =zeros(13, 1); %z-score
        h = zeros(13, 1); %lillietest
       for ii=1:size(list1102,1)
           list=['20181031'; '20181101';'20181102'];
           list1102=['ch09c1';'ch19c1';'ch19c2';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch31c1';'ch43c1';'ch43c2';'ch45c1';'ch45c2']; 
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
                sh=NaN(10000,dim1*dim2);
                sh(:,v)=vmp.data.maps_adsmsh(:,v);
                
                dim1 = 40;
                dim2 = 40;
                smooth = [address,'\smoothed_'];
                %for ise1pixel3
                    if ii>1
                        load(filename)
                    end
                    ise_out = ise1pixel3(actual, sh, dim1, dim2);
                    ise(ii) = ise_out(1,:);
                    ise_sh(:,1,ii) = ise_out(2:end,:);
                   
                    ise_sh1 = ise_sh(:,1,ii);
                    ise_sh1(ise_sh1==0) = []; %exclude zero ise
                    
                    ise_2_5(ii, 1) = prctile(ise_sh1, 2.5);
                    ise_97(ii, 1) = prctile(ise_sh1, 97.5);
                    
                    z(ii, 1) =(ise(ii)-mean(ise_sh1))/std(ise_sh1);
                    
                    h(ii, 1)=lillietest(ise_sh1); %1=skewed. 0=normallly distributed
                   
                    %save new .mat file
                    if ii==1
                        filename = [smooth 'ise1pixel3_p50_15.mat'];
                    end
                    save(filename,'origin','ise','ise_2_5','ise_97','ise_sh','z','h');
            cd ..
       end
       cd .. %back one directory
    end
end
