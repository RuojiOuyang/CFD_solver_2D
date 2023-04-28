clear
clc

data=importdata('naca0012_09_17.txt');%读文件
Jmax=data(3);
Imax=data(2);
num=data(2)*data(3);%节点数量
fid=fopen('naca0012_09_17data.txt','wt');
fprintf(fid,'%d ',Imax);
fprintf(fid,'%d\n',Jmax);

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
        fprintf(fid,'%.16f ',data1(i,j,1) + 0.25);
        fprintf(fid,'%.16f\n',data1(i,j,2));
%         plot(data1(i,j,1) + 0.25,data1(i,j,2),'.r');
%         hold on;
    end
end
fclose(fid);