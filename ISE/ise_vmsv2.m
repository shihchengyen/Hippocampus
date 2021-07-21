% view map
cwd ='C:\Users\Teddy\Downloads\data\New folder';
cd(cwd);
list=['20181031'; '20181101';'20181102'];

for i=1:size(list,1)
    list=['20181031'; '20181101';'20181102'];
    cd(list(i,:));
    if i==1 %1031
       list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
%        start from ch35ca
       for ii=1:size(list1031,1)
            list1031=['ch19c1';'ch19c2';'ch19c3';'ch26c1';'ch26c2';'ch29c1';'ch30c1';'ch30c2';'ch35c1';'ch35c2';'ch35c3';'ch43c1';'ch43c2';'ch45c1';'ch45c2'];
            cd(list1031(ii,:));
            
            %ISE
            clear
            load('vmsv.mat')

            % Canvas padding
            % create overall map and insert padded portions in, to account for
            % cross-portion pairs
            firing_rates = [vms.data.maps_adsm ; vms.data.maps_adsmsh];
            canvas = nan(51, 161, size(vms.data.maps_adsmsh,1) + 1); 
            tic;
            
            % flooring
            floor_padded = nan(42,42,size(vms.data.maps_adsmsh,1)+1);
            floor_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);
            floor_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,3203:3203+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1),1);
            floor_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,3243:3243+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1);
            floor_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,3283:3283+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1);
            floor_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,3323:3323+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1), 2);
            canvas(10:end,1:42,:) = floor_padded;
            
            % ceiling
            ceiling_padded = nan(42,42,size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,1603:3202),size(firing_rates,1),40,40), [3 2 1]), 1);
            ceiling_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,4323:4323+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1),1);
            ceiling_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,4363:4363+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,4403:4403+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,4443:4443+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1), 2);
            canvas(10:end,44:85,:) = ceiling_padded;
                        
            % walls
            walls_padded = nan(8,161,size(vms.data.maps_adsmsh,1)+1);
            walls_padded(:,1:end-1,:) = flip(permute(reshape(firing_rates(:,3203:3203+1280-1), size(vms.data.maps_adsmsh,1)+1, 40*4, 8),[3 2 1]), 1);
            walls_padded(:,end,:) = walls_padded(:,1,:);
            canvas(1:8,:,:) = walls_padded;
            
            % used to pad pillar base more easily
            floor_base = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);

            % pillars
            PTL_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PTL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4963:4963+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            % small diagonal issue here, diagonal floor bins at the corners are put
            % side by side, only 16 such occurrences in total, neglected for now.
            PTL_padded(end,1:8,:) = flip(permute(floor_base(9:16,8,:),[2 1 3]),2);
            PTL_padded(end,9:16,:) = floor_base(8,9:16,:);
            PTL_padded(end,17:24,:) = permute(floor_base(9:16,17,:),[2 1 3]);
            PTL_padded(end,25:32,:) = flip(floor_base(17,9:16,:),2);
            PTL_padded(:,end,:) = PTL_padded(:,1,:);
            canvas(10:10+6-1,87:87+32,:) = PTL_padded;
            
            PTR_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PTR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4803:4803+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PTR_padded(end,1:8,:) = flip(permute(floor_base(9:16,24,:),[2 1 3]),2);
            PTR_padded(end,9:16,:) = floor_base(8,25:32,:);
            PTR_padded(end,17:24,:) = permute(floor_base(9:16,33,:),[2 1 3]);
            PTR_padded(end,25:32,:) = flip(floor_base(17,25:32,:),2);
            PTR_padded(:,end,:) = PTR_padded(:,1,:);
            canvas(10:10+6-1,121:121+32,:) = PTR_padded;
            
            PBL_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PBL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4643:4643+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PBL_padded(end,1:8,:) = flip(permute(floor_base(25:32,8,:),[2 1 3]),2);
            PBL_padded(end,9:16,:) = floor_base(24,9:16,:);
            PBL_padded(end,17:24,:) = permute(floor_base(25:32,17,:),[2 1 3]);
            PBL_padded(end,25:32,:) = flip(floor_base(33,9:16,:),2);
            PBL_padded(:,end,:) = PBL_padded(:,1,:);
            canvas(17:17+6-1,87:87+32,:) = PBL_padded;
            
            PBR_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PBR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4483:4483+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PBR_padded(end,1:8,:) = flip(permute(floor_base(25:32,24,:),[2 1 3]),2);
            PBR_padded(end,9:16,:) = floor_base(24,25:32,:);
            PBR_padded(end,17:24,:) = permute(floor_base(25:32,33,:),[2 1 3]);
            PBR_padded(end,25:32,:) = flip(floor_base(33,25:32,:),2);
            PBR_padded(:,end,:) = PBR_padded(:,1,:);
            canvas(17:17+6-1,121:121+32,:) = PBR_padded;
            
            actual_image = canvas(:,:,1);
            actual_image = actual_image(:)';
            shuffled_images = canvas(:,:,2:end);
            shuffled_images = reshape(shuffled_images, size(shuffled_images,3),size(shuffled_images,1)*size(shuffled_images,2));
            
            disp(['time taken to pad map for ISE: ' num2str(toc)]);
            
            %ISE                          
            %for ise13.mat            
            ise_out1 = ise13(reshape(floor_padded(:,:,1),1,[]), permute(reshape(floor_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 42, 42);
            ise_out2 = ise13(reshape(ceiling_padded(:,:,1),1,[],1), permute(reshape(ceiling_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 42, 42);
            ise_out3 = ise13(reshape(walls_padded(:,:,1),1,[],1), permute(reshape(walls_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 8, 161);
            ise_out4 = ise13(reshape(PTL_padded(:,:,1),1,[],1), permute(reshape(PTL_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out5 = ise13(reshape(PTR_padded(:,:,1),1,[],1), permute(reshape(PTR_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out6 = ise13(reshape(PBL_padded(:,:,1),1,[],1), permute(reshape(PBL_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out7 = ise13(reshape(PBR_padded(:,:,1),1,[],1), permute(reshape(PBR_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out = ise_out1 + ise_out2 + ise_out3 + ise_out4 +ise_out5+ise_out6+ise_out7;      
                
            
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                ise_97 = prctile(ise_sh, 97.5);
                z =(ise-mean(ise_sh))/std(ise_sh); 
                %save new .mat file
                save('ise13_s.mat','ise','ise_2_5','ise_97','z','ise_sh');
                
            
            
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
            clear
            load('vmsv.mat')

            % Canvas padding
            % create overall map and insert padded portions in, to account for
            % cross-portion pairs
            firing_rates = [vms.data.maps_adsm ; vms.data.maps_adsmsh];
            canvas = nan(51, 161, size(vms.data.maps_adsmsh,1) + 1); 
            tic;
            
            % flooring
            floor_padded = nan(42,42,size(vms.data.maps_adsmsh,1)+1);
            floor_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);
            floor_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,3203:3203+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1),1);
            floor_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,3243:3243+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1);
            floor_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,3283:3283+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1);
            floor_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,3323:3323+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1), 2);
            canvas(10:end,1:42,:) = floor_padded;
            
            % ceiling
            ceiling_padded = nan(42,42,size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,1603:3202),size(firing_rates,1),40,40), [3 2 1]), 1);
            ceiling_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,4323:4323+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1),1);
            ceiling_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,4363:4363+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,4403:4403+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,4443:4443+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1), 2);
            canvas(10:end,44:85,:) = ceiling_padded;
                        
            % walls
            walls_padded = nan(8,161,size(vms.data.maps_adsmsh,1)+1);
            walls_padded(:,1:end-1,:) = flip(permute(reshape(firing_rates(:,3203:3203+1280-1), size(vms.data.maps_adsmsh,1)+1, 40*4, 8),[3 2 1]), 1);
            walls_padded(:,end,:) = walls_padded(:,1,:);
            canvas(1:8,:,:) = walls_padded;
            
            % used to pad pillar base more easily
            floor_base = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);

            % pillars
            PTL_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PTL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4963:4963+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            % small diagonal issue here, diagonal floor bins at the corners are put
            % side by side, only 16 such occurrences in total, neglected for now.
            PTL_padded(end,1:8,:) = flip(permute(floor_base(9:16,8,:),[2 1 3]),2);
            PTL_padded(end,9:16,:) = floor_base(8,9:16,:);
            PTL_padded(end,17:24,:) = permute(floor_base(9:16,17,:),[2 1 3]);
            PTL_padded(end,25:32,:) = flip(floor_base(17,9:16,:),2);
            PTL_padded(:,end,:) = PTL_padded(:,1,:);
            canvas(10:10+6-1,87:87+32,:) = PTL_padded;
            
            PTR_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PTR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4803:4803+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PTR_padded(end,1:8,:) = flip(permute(floor_base(9:16,24,:),[2 1 3]),2);
            PTR_padded(end,9:16,:) = floor_base(8,25:32,:);
            PTR_padded(end,17:24,:) = permute(floor_base(9:16,33,:),[2 1 3]);
            PTR_padded(end,25:32,:) = flip(floor_base(17,25:32,:),2);
            PTR_padded(:,end,:) = PTR_padded(:,1,:);
            canvas(10:10+6-1,121:121+32,:) = PTR_padded;
            
            PBL_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PBL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4643:4643+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PBL_padded(end,1:8,:) = flip(permute(floor_base(25:32,8,:),[2 1 3]),2);
            PBL_padded(end,9:16,:) = floor_base(24,9:16,:);
            PBL_padded(end,17:24,:) = permute(floor_base(25:32,17,:),[2 1 3]);
            PBL_padded(end,25:32,:) = flip(floor_base(33,9:16,:),2);
            PBL_padded(:,end,:) = PBL_padded(:,1,:);
            canvas(17:17+6-1,87:87+32,:) = PBL_padded;
            
            PBR_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PBR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4483:4483+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PBR_padded(end,1:8,:) = flip(permute(floor_base(25:32,24,:),[2 1 3]),2);
            PBR_padded(end,9:16,:) = floor_base(24,25:32,:);
            PBR_padded(end,17:24,:) = permute(floor_base(25:32,33,:),[2 1 3]);
            PBR_padded(end,25:32,:) = flip(floor_base(33,25:32,:),2);
            PBR_padded(:,end,:) = PBR_padded(:,1,:);
            canvas(17:17+6-1,121:121+32,:) = PBR_padded;
            
            actual_image = canvas(:,:,1);
            actual_image = actual_image(:)';
            shuffled_images = canvas(:,:,2:end);
            shuffled_images = reshape(shuffled_images, size(shuffled_images,3),size(shuffled_images,1)*size(shuffled_images,2));
            
            disp(['time taken to pad map for ISE: ' num2str(toc)]);
            
            %ISE            
            %for ise13.mat            
            ise_out1 = ise13(reshape(floor_padded(:,:,1),1,[]), permute(reshape(floor_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 42, 42);
            ise_out2 = ise13(reshape(ceiling_padded(:,:,1),1,[],1), permute(reshape(ceiling_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 42, 42);
            ise_out3 = ise13(reshape(walls_padded(:,:,1),1,[],1), permute(reshape(walls_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 8, 161);
            ise_out4 = ise13(reshape(PTL_padded(:,:,1),1,[],1), permute(reshape(PTL_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out5 = ise13(reshape(PTR_padded(:,:,1),1,[],1), permute(reshape(PTR_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out6 = ise13(reshape(PBL_padded(:,:,1),1,[],1), permute(reshape(PBL_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out7 = ise13(reshape(PBR_padded(:,:,1),1,[],1), permute(reshape(PBR_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out = ise_out1 + ise_out2 + ise_out3 + ise_out4 +ise_out5+ise_out6+ise_out7;      
                
            
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                ise_97 = prctile(ise_sh, 97.5);
                z =(ise-mean(ise_sh))/std(ise_sh); 
                %save new .mat file
                save('ise13_s.mat','ise','ise_2_5','ise_97','z','ise_sh');
            
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
            clear
            load('vmsv.mat')

            % Canvas padding
            % create overall map and insert padded portions in, to account for
            % cross-portion pairs
            firing_rates = [vms.data.maps_adsm ; vms.data.maps_adsmsh];
            canvas = nan(51, 161, size(vms.data.maps_adsmsh,1) + 1); 
            tic;
            
            % flooring
            floor_padded = nan(42,42,size(vms.data.maps_adsmsh,1)+1);
            floor_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);
            floor_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,3203:3203+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1),1);
            floor_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,3243:3243+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1);
            floor_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,3283:3283+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1);
            floor_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,3323:3323+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1), 2);
            canvas(10:end,1:42,:) = floor_padded;
            
            % ceiling
            ceiling_padded = nan(42,42,size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,1603:3202),size(firing_rates,1),40,40), [3 2 1]), 1);
            ceiling_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,4323:4323+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1),1);
            ceiling_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,4363:4363+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,4403:4403+39),[2 1]), 40, 1, size(vms.data.maps_adsmsh,1)+1);
            ceiling_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,4443:4443+39),[2 1]), 1, 40, size(vms.data.maps_adsmsh,1)+1), 2);
            canvas(10:end,44:85,:) = ceiling_padded;
                        
            % walls
            walls_padded = nan(8,161,size(vms.data.maps_adsmsh,1)+1);
            walls_padded(:,1:end-1,:) = flip(permute(reshape(firing_rates(:,3203:3203+1280-1), size(vms.data.maps_adsmsh,1)+1, 40*4, 8),[3 2 1]), 1);
            walls_padded(:,end,:) = walls_padded(:,1,:);
            canvas(1:8,:,:) = walls_padded;
            
            % used to pad pillar base more easily
            floor_base = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);

            % pillars
            PTL_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PTL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4963:4963+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            % small diagonal issue here, diagonal floor bins at the corners are put
            % side by side, only 16 such occurrences in total, neglected for now.
            PTL_padded(end,1:8,:) = flip(permute(floor_base(9:16,8,:),[2 1 3]),2);
            PTL_padded(end,9:16,:) = floor_base(8,9:16,:);
            PTL_padded(end,17:24,:) = permute(floor_base(9:16,17,:),[2 1 3]);
            PTL_padded(end,25:32,:) = flip(floor_base(17,9:16,:),2);
            PTL_padded(:,end,:) = PTL_padded(:,1,:);
            canvas(10:10+6-1,87:87+32,:) = PTL_padded;
            
            PTR_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PTR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4803:4803+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PTR_padded(end,1:8,:) = flip(permute(floor_base(9:16,24,:),[2 1 3]),2);
            PTR_padded(end,9:16,:) = floor_base(8,25:32,:);
            PTR_padded(end,17:24,:) = permute(floor_base(9:16,33,:),[2 1 3]);
            PTR_padded(end,25:32,:) = flip(floor_base(17,25:32,:),2);
            PTR_padded(:,end,:) = PTR_padded(:,1,:);
            canvas(10:10+6-1,121:121+32,:) = PTR_padded;
            
            PBL_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PBL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4643:4643+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PBL_padded(end,1:8,:) = flip(permute(floor_base(25:32,8,:),[2 1 3]),2);
            PBL_padded(end,9:16,:) = floor_base(24,9:16,:);
            PBL_padded(end,17:24,:) = permute(floor_base(25:32,17,:),[2 1 3]);
            PBL_padded(end,25:32,:) = flip(floor_base(33,9:16,:),2);
            PBL_padded(:,end,:) = PBL_padded(:,1,:);
            canvas(17:17+6-1,87:87+32,:) = PBL_padded;
            
            PBR_padded = nan(6,33,size(vms.data.maps_adsmsh,1)+1);
            PBR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4483:4483+160-1), size(vms.data.maps_adsmsh,1)+1, 8*4, 5),[3 2 1]), 1);
            PBR_padded(end,1:8,:) = flip(permute(floor_base(25:32,24,:),[2 1 3]),2);
            PBR_padded(end,9:16,:) = floor_base(24,25:32,:);
            PBR_padded(end,17:24,:) = permute(floor_base(25:32,33,:),[2 1 3]);
            PBR_padded(end,25:32,:) = flip(floor_base(33,25:32,:),2);
            PBR_padded(:,end,:) = PBR_padded(:,1,:);
            canvas(17:17+6-1,121:121+32,:) = PBR_padded;
            
            actual_image = canvas(:,:,1);
            actual_image = actual_image(:)';
            shuffled_images = canvas(:,:,2:end);
            shuffled_images = reshape(shuffled_images, size(shuffled_images,3),size(shuffled_images,1)*size(shuffled_images,2));
            
            disp(['time taken to pad map for ISE: ' num2str(toc)]);
            
            %ISE            
            %for ise13.mat            
            ise_out1 = ise13(reshape(floor_padded(:,:,1),1,[]), permute(reshape(floor_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 42, 42);
            ise_out2 = ise13(reshape(ceiling_padded(:,:,1),1,[],1), permute(reshape(ceiling_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 42, 42);
            ise_out3 = ise13(reshape(walls_padded(:,:,1),1,[],1), permute(reshape(walls_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 8, 161);
            ise_out4 = ise13(reshape(PTL_padded(:,:,1),1,[],1), permute(reshape(PTL_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out5 = ise13(reshape(PTR_padded(:,:,1),1,[],1), permute(reshape(PTR_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out6 = ise13(reshape(PBL_padded(:,:,1),1,[],1), permute(reshape(PBL_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out7 = ise13(reshape(PBR_padded(:,:,1),1,[],1), permute(reshape(PBR_padded(:,:,2:end),1,[],size(vms.data.maps_adsmsh,1)),[3 2 1]), 6, 33);
            ise_out = ise_out1 + ise_out2 + ise_out3 + ise_out4 +ise_out5+ise_out6+ise_out7;      
                
             
                ise = ise_out(1);
                ise_sh = ise_out(2:end);
                ise_sh(ise_sh==0) = []; %exclude zero ise
                ise_2_5 = prctile(ise_sh, 2.5);
                ise_97 = prctile(ise_sh, 97.5);
                z =(ise-mean(ise_sh))/std(ise_sh); 
                %save new .mat file
                save('ise13_s.mat','ise','ise_2_5','ise_97','z','ise_sh');   
           
            
            
            cd ..
       end
       cd .. %back one directory
    end
end

