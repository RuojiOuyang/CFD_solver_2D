data=importdata('yuanzhu.txt');%读文件
Jmax=data(2);
Imax=data(3);
num=data(2)*data(3);%节点数量
fid=fopen('yuanzhudata.txt','wt');
fprintf(fid,'%d ',data(3));
fprintf(fid,'%d\n',data(2));
k=4;%这里的4是由网格数据决定的
for i=1:Imax
    for j=1:Jmax
        data1(i,j,1)=data(k);
        data1(i,j,2)=data(k+num);
        k=k+1;
    end
end

for i=1:Imax
    for j=1:Jmax
        data2(i,j,1)=data1(Imax-i+1,j,1);
        data2(i,j,2)=data1(Imax-i+1,j,2);
    end
end

for j=1:Jmax
    for i=1:Imax
        fprintf(fid,'%.16f ',data2(i,j,1));
        fprintf(fid,'%.16f\n',data2(i,j,2));
%         plot(data2(i,j,1),data2(i,j,2),'.r');
%         hold on;
    end
end
fclose(fid);