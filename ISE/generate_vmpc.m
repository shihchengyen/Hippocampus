clear
clc

cwd ='C:\Users\Teddy\Downloads\data'; %ASUS
cd(cwd);
list=['20181031\session01\array01'; '20181101\session01\array01';'20181102\session01\array01'];

for i=1:size(list,1)
    list=['20181031\session01\array01'; '20181101\session01\array01';'20181102\session01\array01'];
    
    cd(list(i,:));
    if i==1 %1031
       list1031=['channel019\cell01';'channel019\cell02';'channel019\cell03';'channel026\cell01';'channel026\cell02';'channel029\cell01';'channel030\cell01';'channel030\cell02';'channel035\cell01';'channel035\cell02';'channel035\cell03';'channel043\cell01';'channel043\cell02';'channel045\cell01';'channel045\cell02'];
       
       for ii=1:size(list1031,1)
            list1031=['channel019\cell01';'channel019\cell02';'channel019\cell03';'channel026\cell01';'channel026\cell02';'channel029\cell01';'channel030\cell01';'channel030\cell02';'channel035\cell01';'channel035\cell02';'channel035\cell03';'channel043\cell01';'channel043\cell02';'channel045\cell01';'channel045\cell02'];
            cd(list1031(ii,:));
%             vmp=vmpc('auto');save('vmpc_ted','vmp');
            vms=vmsv2('auto');save('vmsv_ted','vms');
            
            cd ..
            cd ..
       end
       cd .. %back one directory
    end
    
    if i==2 %1101
       list1101=['channel019\cell01';'channel019\cell02';'channel021\cell01';'channel023\cell01';'channel029\cell01';'channel029\cell02';'channel029\cell03';'channel029\cell04';'channel030\cell01';'channel030\cell02';'channel035\cell01';'channel043\cell01';'channel045\cell01'];
       for ii=1:size(list1101,1)
           list1101=['channel019\cell01';'channel019\cell02';'channel021\cell01';'channel023\cell01';'channel029\cell01';'channel029\cell02';'channel029\cell03';'channel029\cell04';'channel030\cell01';'channel030\cell02';'channel035\cell01';'channel043\cell01';'channel045\cell01']; 
           cd(list1101(ii,:));
%            vmp=vmpc('auto');save('vmpc_ted','vmp');
           vms=vmsv2('auto');save('vmsv_ted','vms');
           cd ..
           cd ..
       end
       cd .. %back one directory
    end
    
    if i==3 %1102
       list1102=['channel009\cell01';'channel019\cell01';'channel019\cell02';'channel026\cell01';'channel026\cell02';'channel029\cell01';'channel030\cell01';'channel030\cell02';'channel031\cell01';'channel043\cell01';'channel043\cell02';'channel045\cell01';'channel045\cell02'];
       for ii=1:size(list1102,1)
           list1102=['channel009\cell01';'channel019\cell01';'channel019\cell02';'channel026\cell01';'channel026\cell02';'channel029\cell01';'channel030\cell01';'channel030\cell02';'channel031\cell01';'channel043\cell01';'channel043\cell02';'channel045\cell01';'channel045\cell02']; 
           cd(list1102(ii,:));
%            vmp=vmpc('auto');save('vmpc_ted','vmp');
           vms=vmsv2('auto');save('vmsv_ted','vms');
           cd ..
           cd ..
       end
       cd .. %back one directory
    end
end
