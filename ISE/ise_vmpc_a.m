% collect the ise_vmpc
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
                
                %for ise1_1.mat
                load('ise1_1.mat');
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z'); 

                %for ise1_2.mat
                load('ise1_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise1_3.mat
                load('ise1_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise2_1.mat
                load('ise2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise2_2.mat
                load('ise2_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise2_3.mat
                load('ise2_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_1.mat
                load('ise3_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_2.mat
                load('ise3_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3.mat
                load('ise3_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_2_1.mat
                load('ise3_2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_2_2.mat
                load('ise3_2_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_2_3.mat
                load('ise3_2_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3_1.mat
                load('ise3_3_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_3_2.mat
                load('ise3_3_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3_3.mat
                load('ise3_3_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise4_1.mat
                load('ise4_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise4_2.mat
                load('ise4_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise4_3.mat
                load('ise4_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_1.mat
                load('ise5_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_2.mat
                load('ise5_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_3.mat
                load('ise5_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise6_1.mat
                load('ise6_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise6_2.mat
                load('ise6_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise6_3.mat
                load('ise6_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise6_2_1.mat
                load('ise6_2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
            
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
                
                %for ise1_1.mat
                load('ise1_1.mat');
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z'); 

                %for ise1_2.mat
                load('ise1_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise1_3.mat
                load('ise1_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise2_1.mat
                load('ise2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise2_2.mat
                load('ise2_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise2_3.mat
                load('ise2_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_1.mat
                load('ise3_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_2.mat
                load('ise3_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3.mat
                load('ise3_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_2_1.mat
                load('ise3_2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_2_2.mat
                load('ise3_2_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_2_3.mat
                load('ise3_2_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3_1.mat
                load('ise3_3_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_3_2.mat
                load('ise3_3_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3_3.mat
                load('ise3_3_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise4_1.mat
                load('ise4_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise4_2.mat
                load('ise4_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise4_3.mat
                load('ise4_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_1.mat
                load('ise5_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_2.mat
                load('ise5_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_3.mat
                load('ise5_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise6_1.mat
                load('ise6_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise6_2.mat
                load('ise6_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise6_3.mat
                load('ise6_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise6_2_1.mat
                load('ise6_2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
            
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
                
                %for ise1_1.mat
                load('ise1_1.mat');
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z'); 

                %for ise1_2.mat
                load('ise1_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise1_3.mat
                load('ise1_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise2_1.mat
                load('ise2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise2_2.mat
                load('ise2_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise2_3.mat
                load('ise2_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_1.mat
                load('ise3_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_2.mat
                load('ise3_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3.mat
                load('ise3_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_2_1.mat
                load('ise3_2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_2_2.mat
                load('ise3_2_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_2_3.mat
                load('ise3_2_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3_1.mat
                load('ise3_3_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise3_3_2.mat
                load('ise3_3_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise3_3_3.mat
                load('ise3_3_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise4_1.mat
                load('ise4_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise4_2.mat
                load('ise4_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise4_3.mat
                load('ise4_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_1.mat
                load('ise5_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_2.mat
                load('ise5_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise5_3.mat
                load('ise5_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise6_1.mat
                load('ise6_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise6_2.mat
                load('ise6_2.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 

                %for ise6_3.mat
                load('ise6_3.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append'); 
                
                %for ise6_2_1.mat
                load('ise6_2_1.mat');
                Origin = origin;
                ISE = ise;
                ISE_2_5 = ise_2_5;
                Z =z;
                load('ise_p.mat');
                origin = [origin;Origin]; 
                ise = [ise; ISE];
                ise_2_5 = [ise_2_5;ISE_2_5];
                z = [z;Z];
                %save new .mat file
                save('ise_p.mat','origin','ise','ise_2_5','z','-append');    
             
            
            cd ..
       end
       cd .. %back one directory
    end
end

