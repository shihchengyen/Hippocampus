function t = hyclusterquality(thresholdi,thresholdl,thresholdr)
cd('output')
moutput = readmda('firings.mda');
cd('..')
cd('dataset')
morigin = readmda('raw_data.mda');
cd('..')
points=[];
row0=1;
column0=1;
row1=1;
column1=1;
maxcluster=max(moutput(3,:));
isodistance=[];
lratio=[];
interval=cell(1,maxcluster);
rperiod=[];
index=[];
index1=[];
size1=0;
delete=[];
clusterquality=[];
rthreshold=2.5;
tt=cell(1,maxcluster);
twidth=50;
maxx=-100;
swidth=cell(1,maxcluster);
meanwidth=[];
amp=cell(1,maxcluster);
meanamp=[];
widthvar=[];
ampvar=[];
good=0;
if nargin<3
    thresholdi=16;
    thresholdl=0.1;
    thresholdr=0.01;
end
for k = 1:maxcluster
    clusterquality(k,1)=k;
end
for i=moutput(2,:)
    column0=1;
    for j=(i-25):(i+25)
        points(row0,column0)=morigin(j);
        column0=column0+1;
    end
    
        
    row0=row0+1;
end

[coeff, score] =pca(points);
score1=score(:,1:5);
for i = 1:maxcluster
    column1=1;
    for j=1:size(moutput,2)
        if moutput(3,j)==i
            index(row1,column1)=j;
        end
        column1=column1+1;
    end
    row1=row1+1;
end
  for ii =1:maxcluster
    delete=[];
    index1=index(ii,:);
    size1=size(index1,2);
    for i=1:size1
        if index1(i)==0
            delete=[delete,i];
        end
    end
    index1(delete)=[];
    for iii=2:size(index1,2)
        intervaln=moutput(2,index1(iii))-moutput(2,index1(iii-1));
        intervalt=intervaln/30;
        interval{ii}(iii-1)=intervalt;
    end
    IsolDist = IsolationDistance(score1, index1);
    isodistance=[isodistance,IsolDist];
    [L, Lratio, df] = L_Ratio(score1, index1);
    lratio=[lratio,Lratio];
    
  end
isodistance=isodistance';
lratio=lratio';
for i = 1:maxcluster
    rperiodn=sum(interval{1,i}<2.5)/(length(interval{1,i})+1);
    rperiod=[rperiod,rperiodn];
end
rperiod=rperiod';
clusterquality(:,2)=isodistance;
clusterquality(:,3)=lratio;
clusterquality(:,4)=rperiod;
for i =1:maxcluster
    for k = 1:size(moutput,2)
        if moutput(3,k)==i
            tt{1,i}=[tt{1,i},moutput(2,k)];
        end
    end
end
for j=1:maxcluster
    for k=1:size(tt{1,j},2)
        maxx=-100;
        minn=100;
        for i=tt{1,j}(k):tt{1,j}(k)+twidth
            if morigin(tt{j}(k))<0
                if morigin(i)>maxx
                    maxx=morigin(i);
                    tvalue=i;
                end
            else
                if morigin(i)<minn
                    minn=morigin(i);
                    tvalue=i;
                end
            end
        end
        swidth{j}=[swidth{j},tvalue-tt{j}(k)];
    end
end
for i = 1:maxcluster
    meanwidth(i)=mean(swidth{i});
end
meanwidth=meanwidth';
clusterquality(:,5)=meanwidth;
for i =1:maxcluster
    amp{i}=morigin(tt{i});
    meanamp=[meanamp,mean(amp{i})];
end
meanamp=meanamp';
clusterquality(:,6)=meanamp;
for i =1:maxcluster
    widthvar(i,1)=var(swidth{i});
    ampvar(i,1)=var(amp{i});
end
clusterquality(:,7)=widthvar;
clusterquality(:,8)=ampvar;
for i = 1:maxcluster
    good=0;
    if clusterquality(i,2)>thresholdi
        good=good+1;
    end
    if clusterquality(i,3)<thresholdl
        good=good+1;
    end
    if clusterquality(i,4)<thresholdr
        good=good+1;
    end
    clusterquality(i,9)=good;
end
overallscore=clusterquality(:,9);
clusterindex=clusterquality(:,1);
t=table(clusterindex,isodistance,lratio,rperiod,meanwidth,meanamp,widthvar,ampvar,overallscore);
save('clusterquality.mat','t')
isodistancesf=round(isodistance,3,'significant');
lratiosf=round(lratio,3,'significant');
rperiodsf=round(rperiod,3,'significant');
meanwidthsf=round(meanwidth,3,'significant');
meanampsf=round(meanamp,3,'significant');
widthvarsf=round(widthvar,3,'significant');
ampvarsf=round(ampvar,3,'significant');
tsf=table(clusterindex,isodistancesf,lratiosf,rperiodsf,meanwidthsf,meanampsf,widthvarsf,ampvarsf,overallscore);
writetable(tsf,'clusterquality.csv')

end