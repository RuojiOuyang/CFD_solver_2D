clear
clc

data=importdata('naca0012.txt');%读文件
Jmax=data(2);
Imax=data(3);
num=data(2)*data(3);%节点数量
fid=fopen('naca0012data.txt','wt');
fprintf(fid,'%d ',data(3));
fprintf(fid,'%d\n',data(2));

data1 = zeros(Imax,Jmax,2);

k=5;%这里的4是由网格数据决定的
for j=1:Jmax
    for i=1:Imax
        data1(i,j,1)=data(k);
        data1(i,j,2)=data(k+num);
        k=k+1;
    end
end

for j=1:Jmax
    for i=1:Imax
        fprintf(fid,'%.16f ',data1(i,j,1));
        fprintf(fid,'%.16f\n',data1(i,j,2));
%         plot(data1(i,j,1),data1(i,j,2),'.r');
%         hold on;
    end
end
fclose(fid);