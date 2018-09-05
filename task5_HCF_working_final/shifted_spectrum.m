%shift correction
close all
clear all
n0=1;
I_const=0.5*n0*const.c*const.eps0;


Inpt_range=[1;1001];%[1175;1683];
Outpt_range=[1;1001];
[x_data,y_data,x_data2,y_data2]=extract_data('shiftMe.fig',Inpt_range,Outpt_range);%1175:1683


shifted_ydata=zeros(1,length(y_data));
shifted_ydata(1:989)=y_data(13:1001);

[lambda,shifted_IntensityL]=convert_data(x_data,shifted_ydata,'back');
[lambda2,hcf_IntensityL]=convert_data(x_data,y_data2,'back');

[lambda,shifted_IntensityL]=even_data(lambda,shifted_IntensityL);
[lambda2,hcf_IntensityL]=even_data(lambda2,hcf_IntensityL);

figure;
subplot(2,1,1)
plot(x_data,[y_data2;shifted_ydata])
hold on
subplot(2,1,2)
%plot(lambda,[hcf_IntensityL;shifted_IntensityL])
plot(lambda(1:100).*1e9,[hcf_IntensityL(1:100);shifted_IntensityL(1:100)])