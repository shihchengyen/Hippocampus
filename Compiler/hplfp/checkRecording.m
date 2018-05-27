function checkRecording

hd = ProcessLevel(hidata,'Levels','Session','FileName','rplhighpass.mat');

figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,1,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'ha1.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,2,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'ha2.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,3,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'ha3.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,4,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'ha4.png')
close

%{
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,5,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs2a1.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,6,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs2a2.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,7,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs2a3.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,8,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs2a4.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,9,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs3a1.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,10,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs3a2.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,11,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs3a3.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,12,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs3a4.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,13,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs4a1.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,14,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs4a2.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,15,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs4a3.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,16,'Array','UseGMR','Objects',{'rplhighpass',{'FFT'},{'HighpassFreqs',[500 10000]}});
saveas(gcf,'hs4a4.png')
close
%}

figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,1,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'la1.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,2,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'la2.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,3,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'la3.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,4,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'la4.png')
close

%{
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,5,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls2a1.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,6,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls2a2.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,7,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls2a3.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,8,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls2a4.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,9,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls3a1.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,10,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls3a2.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,11,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls3a3.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,12,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls3a4.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,13,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls4a1.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,14,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls4a2.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,15,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls4a3.png')
close
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(hd,16,'Array','UseGMR','Objects',{'rpllfp',{'FFT','XLims',[0 150]}});
saveas(gcf,'ls4a4.png')
close
%}
