data1=importdata('result.dat');

Imax=197;
Jmax=30;

num=1;
for j=1:1:Jmax
    for i=1:1:Imax
        x(i,j)=data1.data(num,1);
        y(i,j)=data1.data(num,2);
        Den(i,j)=data1.data(num,3);
        p(i,j)=data1.data(num,4);
        T(i,j)=data1.data(num,5);
        Ma(i,j)=data1.data(num,6);
        u(i,j)=data1.data(num,7);
        v(i,j)=data1.data(num,8);
        u_y(i,j)=data1.data(num,11);
        v_x(i,j)=data1.data(num,12);
        num=num+1;
    end
%     x(Imax+1,j)=x(1,j);
%     y(Imax+1,j)=y(1,j);
%     Den(Imax+1,j)=Den(1,j);
%     p(Imax+1,j)=p(1,j);
%     T(Imax+1,j)=T(1,j);
%     Ma(Imax+1,j)=Ma(1,j);
end
wz=v_x-u_y;
wz=wz.*wz;
figure('name','ÃÜ¶ÈÔÆÍ¼')
pcolor(x,y,Den);
colormap jet,shading interp;
figure('name','Ñ¹Á¦ÔÆÍ¼')
pcolor(x,y,p);
colormap jet,shading interp;
figure('name','ÎÂ¶ÈÔÆÍ¼')
pcolor(x,y,T);
colormap jet,shading interp;
figure('name','ÂíºÕÊıÔÆÍ¼')
pcolor(x,y,Ma);
colormap jet,shading interp;
figure('name','ÎĞÁ¿ÔÆÍ¼')
pcolor(x,y,wz);
colormap jet,shading interp;