data1=importdata('SteadyInitialization0.dat');

Imax=181;
Jmax=90;

num=1;
for j=1:1:Jmax
    for i=1:1:Imax
        x(i,j)=data1(num,1);
        y(i,j)=data1(num,2);
        Den(i,j)=data1(num,3);
        p(i,j)=data1(num,4);
        T(i,j)=data1(num,5);
        Ma(i,j)=data1(num,6);
        num=num+1;
    end
    x(Imax+1,j)=x(1,j);
    y(Imax+1,j)=y(1,j);
    Den(Imax+1,j)=Den(1,j);
    p(Imax+1,j)=p(1,j);
    T(Imax+1,j)=T(1,j);
    Ma(Imax+1,j)=Ma(1,j);
end

figure('name','等压线图')
MaxP=floor(max(max(p)));
MinP=ceil(min(min(p)));
delta=floor((MaxP-MinP)/10);
WriteP=MinP:delta:MaxP;
[cs, h]=contour(x,y,p,WriteP,'k','LineWidth',1);
axis([-0.5 1.5 -1 1]);
clabel(cs,h,'FontSize',10,'Color','k','Rotation',0,'LabelSpacing',1000);
hold on;

x1=0:0.01:1;
y1=0.6*(-0.1015*x1.^4+0.2843*x1.^3-0.3576*x1.^2-0.1221*x1+0.2969*sqrt(x1));
y2=-y1;
plot(x1,y1,'k','LineWidth',1);
hold on;
plot(x1,y2,'k','LineWidth',1);

data = importdata('Res.dat');%读文件
figure('name','残差')
plot(data(:,1),data(:,2),'k');