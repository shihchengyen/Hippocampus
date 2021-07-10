cwd ='C:\Users\Teddy\Downloads\data\New folder';
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
                
                %for ise1_1.mat
                origin = vmp.data.origin;
                ise_out = ise1_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_1.mat','origin','ise','ise_2_5','z'); 

                %for ise1_2.mat
                origin = vmp.data.origin;
                ise_out = ise1_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise1_3.mat
                origin = vmp.data.origin;
                ise_out = ise1_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_3.mat','origin','ise','ise_2_5','z'); 

                %for ise2_1.mat
                origin = vmp.data.origin;
                ise_out = ise2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise2_2.mat
                origin = vmp.data.origin;
                ise_out = ise2_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_2.mat','origin','ise','ise_2_5','z'); 

                %for ise2_3.mat
                origin = vmp.data.origin;
                ise_out = ise2_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_1.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_2_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_2.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_1.mat','origin','ise','ise_2_5','z'); 

                %for ise3_3_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_3.mat','origin','ise','ise_2_5','z'); 

                %for ise4_1.mat
                origin = vmp.data.origin;
                ise_out = ise4_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise4_2.mat
                origin = vmp.data.origin;
                ise_out = ise4_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_2.mat','origin','ise','ise_2_5','z'); 

                %for ise4_3.mat
                origin = vmp.data.origin;
                ise_out = ise4_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise5_1.mat
                origin = vmp.data.origin;
                ise_out = ise5_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_1.mat','origin','ise','ise_2_5','z'); 

                %for ise5_2.mat
                origin = vmp.data.origin;
                ise_out = ise5_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise5_3.mat
                origin = vmp.data.origin;
                ise_out = ise5_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_3.mat','origin','ise','ise_2_5','z'); 

                %for ise6_1.mat
                origin = vmp.data.origin;
                ise_out = ise6_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise6_2.mat
                origin = vmp.data.origin;
                ise_out = ise6_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_2.mat','origin','ise','ise_2_5','z'); 

                %for ise6_3.mat
                origin = vmp.data.origin;
                ise_out = ise6_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise6_2_1.mat
                origin = vmp.data.origin;
                ise_out = ise6_2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_2_1.mat','origin','ise','ise_2_5','z'); 
            
            cd ..
       end
       cd .. %back one directory
    end
    
    if i==2 %1101
       list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1'];
       for ii=1:size(list1101,1)
           list1101=['ch19c1';'ch19c2';'ch21c1';'ch23c1';'ch29c1';'ch29c2';'ch29c3';'ch30c1';'ch30c2';'ch35c1';'ch43c1';'ch45c1']; 
           cd(list1101(ii,:));
            
            %ISE
                
                clear all
                load('vmpc.mat');   
                
                %for ise1_1.mat
                origin = vmp.data.origin;
                ise_out = ise1_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_1.mat','origin','ise','ise_2_5','z'); 

                %for ise1_2.mat
                origin = vmp.data.origin;
                ise_out = ise1_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise1_3.mat
                origin = vmp.data.origin;
                ise_out = ise1_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_3.mat','origin','ise','ise_2_5','z'); 

                %for ise2_1.mat
                origin = vmp.data.origin;
                ise_out = ise2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise2_2.mat
                origin = vmp.data.origin;
                ise_out = ise2_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_2.mat','origin','ise','ise_2_5','z'); 

                %for ise2_3.mat
                origin = vmp.data.origin;
                ise_out = ise2_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_1.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_2_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_2.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_1.mat','origin','ise','ise_2_5','z'); 

                %for ise3_3_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_3.mat','origin','ise','ise_2_5','z'); 

                %for ise4_1.mat
                origin = vmp.data.origin;
                ise_out = ise4_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise4_2.mat
                origin = vmp.data.origin;
                ise_out = ise4_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_2.mat','origin','ise','ise_2_5','z'); 

                %for ise4_3.mat
                origin = vmp.data.origin;
                ise_out = ise4_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise5_1.mat
                origin = vmp.data.origin;
                ise_out = ise5_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_1.mat','origin','ise','ise_2_5','z'); 

                %for ise5_2.mat
                origin = vmp.data.origin;
                ise_out = ise5_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise5_3.mat
                origin = vmp.data.origin;
                ise_out = ise5_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_3.mat','origin','ise','ise_2_5','z'); 

                %for ise6_1.mat
                origin = vmp.data.origin;
                ise_out = ise6_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise6_2.mat
                origin = vmp.data.origin;
                ise_out = ise6_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_2.mat','origin','ise','ise_2_5','z'); 

                %for ise6_3.mat
                origin = vmp.data.origin;
                ise_out = ise6_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise6_2_1.mat
                origin = vmp.data.origin;
                ise_out = ise6_2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_2_1.mat','origin','ise','ise_2_5','z'); 
            
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
                
                %for ise1_1.mat
                origin = vmp.data.origin;
                ise_out = ise1_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_1.mat','origin','ise','ise_2_5','z'); 

                %for ise1_2.mat
                origin = vmp.data.origin;
                ise_out = ise1_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise1_3.mat
                origin = vmp.data.origin;
                ise_out = ise1_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise1_3.mat','origin','ise','ise_2_5','z'); 

                %for ise2_1.mat
                origin = vmp.data.origin;
                ise_out = ise2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise2_2.mat
                origin = vmp.data.origin;
                ise_out = ise2_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_2.mat','origin','ise','ise_2_5','z'); 

                %for ise2_3.mat
                origin = vmp.data.origin;
                ise_out = ise2_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise2_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_1.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_2_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_2.mat','origin','ise','ise_2_5','z'); 

                %for ise3_2_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_2_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_2_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3_1.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_1.mat','origin','ise','ise_2_5','z'); 

                %for ise3_3_2.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise3_3_3.mat
                origin = vmp.data.origin;
                ise_out = ise3_3_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise3_3_3.mat','origin','ise','ise_2_5','z'); 

                %for ise4_1.mat
                origin = vmp.data.origin;
                ise_out = ise4_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise4_2.mat
                origin = vmp.data.origin;
                ise_out = ise4_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_2.mat','origin','ise','ise_2_5','z'); 

                %for ise4_3.mat
                origin = vmp.data.origin;
                ise_out = ise4_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise4_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise5_1.mat
                origin = vmp.data.origin;
                ise_out = ise5_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_1.mat','origin','ise','ise_2_5','z'); 

                %for ise5_2.mat
                origin = vmp.data.origin;
                ise_out = ise5_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_2.mat','origin','ise','ise_2_5','z'); 
                
                %for ise5_3.mat
                origin = vmp.data.origin;
                ise_out = ise5_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise5_3.mat','origin','ise','ise_2_5','z'); 

                %for ise6_1.mat
                origin = vmp.data.origin;
                ise_out = ise6_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_1.mat','origin','ise','ise_2_5','z'); 
                
                %for ise6_2.mat
                origin = vmp.data.origin;
                ise_out = ise6_2(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_2.mat','origin','ise','ise_2_5','z'); 

                %for ise6_3.mat
                origin = vmp.data.origin;
                ise_out = ise6_3(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_3.mat','origin','ise','ise_2_5','z'); 
                
                %for ise6_2_1.mat
                origin = vmp.data.origin;
                ise_out = ise6_2_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                z =(ise-mean(ise_sh))/std(ise_sh);
                %save new .mat file
                save('ise6_2_1.mat','origin','ise','ise_2_5','z'); 
            
            cd ..
       end
       cd .. %back one directory
    end
end

% %append ise_1_1.mat
% load('ise_1_1.mat')
% Origin = vmp.data.origin;
% origin=[origin;Origin]; %file directory
% ise_out = ise1_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40); %ISE calculation
% ISE = ise_out(1); 
% ise = [ise; ISE];
% ise_sh = ise_out(2:end); %ISE shuffle
% ise_sh(ise_sh==0) = []; %exclude zero ise
% ISE_2_5 = prctile(ise_sh, 2.5); %ISE 2.5 percentile
% ise_2_5 = [ise_2_5;ISE_2_5];
% Z =(ise-mean(ise_sh))/std(ise_sh); %ISE z score from 
% z = [z;Z];
% save('ise_1_1.mat', 'origin', 'ise','ise_2_5','z','-append'); 
