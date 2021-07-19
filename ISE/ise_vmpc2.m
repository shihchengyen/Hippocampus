cwd ='C:\Users\Teddy\Downloads\data\New folder';
% cwd = 'C:\Users\teddy\Downloads\Data\New folder';
cd(cwd);
list=['20181031'; '20181101';'20181102'];


for i=1:size(list,1)
    list=['20181031'; '20181101';'20181102'];
    
    cd(list(i,:));
    if i==1 %1031
       list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
       for ii=1:size(list1031,1)
            list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
            cd(list1031(ii,:));
            
            %ISE
                
                clear all
                load('vmpc.mat');   
                
                %for ise11.mat
                origin = vmp.data.origin;
                ise_out = ise11(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                ise_97 = prctile(ise_sh, 97.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise11_p.mat','origin','ise','ise_2_5','ise_97','z','ise_sh'); 

                
            
            cd ..
       end
       cd .. %back one directory
    end
    
    if i==2 %1101
       list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch29c4';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1'];
       for ii=1:size(list1101,1)
           list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch29c4';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1']; 
           cd(list1101(ii,:));
            
            %ISE
                
                clear all
                load('vmpc.mat');   
                
                %for ise11.mat
                origin = vmp.data.origin;
                ise_out = ise11(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise11_p.mat','origin','ise','ise_2_5','z','ise_sh'); 

            
            cd ..
       end
       cd .. %back one directory
    end
    
    if i==3 %1102
       list1102=['ch09c1';'ch19c1';'ch19c2';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch31c1';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
       for ii=1:size(list1102,1)
           list1102=['ch09c1';'ch19c1';'ch19c2';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch31c1';'ch43c1';'ch43c2';'ch45c1';'ch45c2']; 
           cd(list1102(ii,:));
            
            %ISE
                
                clear all
                load('vmpc.mat');   
                
                %for ise11.mat
                origin = vmp.data.origin;
                ise_out = ise11(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise11_p.mat','origin','ise','ise_2_5','z','ise_sh'); 

            
            cd ..
       end
       cd .. %back one directory
    end
end
