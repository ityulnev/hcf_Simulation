close all
clear all

% close all
%% Parameters

%Laser
wavelength=800e-9;%[m]
f0=const.c/wavelength;
w0=2*pi*f0;
Q_In=2.2e-3;%[J]
t_fwhm=35e-15;%[s]
tau=t_fwhm/(2*sqrt(log(2)));
Power=Q_In/t_fwhm;%[W]


% Hollow Core Fiber
n0=calculate_refr_index(wavelength);
r_hcf=200e-6;%[m]
area_hcf=pi*(r_hcf)^2;%[m^2]
kfreq=n0*w0/const.c;%1/m
Fluence=Q_In/area_hcf;
alpha=0.7885;


n2=1.4e-24;%5.9/////3.6e-19*1e-4;%[m^2/W]  -> SPM
Beta2=-5*(1e-30);%[s^2/m] -> GVD
%best Overlap @ 2.1 and -4.4


%% Meshing
Lz=100e-2;%m
dz=1e-3;
z=0:dz:Lz;

%% Input Signal from Data
Inpt_range=[1;1900];%[1175;1683];
Outpt_range=[1;1900];
[x_data,y_data,x_data2,y_data2]=extract_data('HCFvsFund.fig',Inpt_range,Outpt_range);%1175:1683
%Into SI-Units [nm]->[m]
x_data=x_data.*1e-9;
x_data2=x_data2.*1e-9;

%% Calculate units based in Real center wavelength
weighted_middle=sum(x_data.*y_data)/sum(y_data);
wavelength=weighted_middle;
n0=calculate_refr_index(wavelength);
f0=const.c/wavelength;
w0=2*pi*f0;

%% Signal Amplitude
I_const=0.5*n0*const.c*const.eps0;
A0_inT=sqrt((Fluence)*(1/(sqrt(pi)*I_const*tau)));

%% Prepare data for Fourier Split Step method
[x_dataF,y_dataF]=process_data(x_data,y_data,I_const);

%cut out data
cutrange=48:1048;%48:1048;%18:1100;%Middle for 800nm is 309:775  ****middle @Index 547-548 ***
f=x_dataF(cutrange);
y_dataF=y_dataF(cutrange);
y_dataF=smooth_sides(y_dataF,460,555);


df=f(2)-f(1);
t=linspace(-1/(2*df),1/(2*df),length(f));%(4*pi).*

scaling=sqrt(Fluence/(sum(get_Intensity(y_dataF,I_const))*df));
signalInF=y_dataF.*scaling;
signalInT=fftshift(ifft(ifftshift(signalInF))).*(length(f)*df);%./(sqrt(4*pi))
% signalInW=fftshift(fft(ifftshift(signalInT)))./(length(w)*dw);
% Lu_best_scale=sum(abs(signalInW).^2.*(0.5*const.c.*const.eps0))*dw/Fluence;
% signalInW=signalInW./sqrt(Lu_best_scale);
hcf_IntensityT=I_const.*abs(signalInT).^2;
hcf_IntensityF=I_const.*abs(signalInF).^2;
% trapz(t,hcf_IntensityT)
% trapz(f,hcf_IntensityF)
% plot(t,hcf_IntensityT)




%% Real Data after HCF propagation
[x_dataF2,y_dataF2]=process_data(x_data2,y_data2,I_const);
y_dataF2=y_dataF2(cutrange);


scaling_Eloss=1/2.2;
scaling2=sqrt(scaling_Eloss*Fluence/(sum(get_Intensity(y_dataF2,I_const))*df));

signalInF2=y_dataF2.*scaling2;
hcf_IntensityF2=I_const.*abs(signalInF2).^2;
signalInT2=fftshift(ifft(ifftshift(signalInF2))).*(length(f)*df);%./sqrt(4*pi)
hcf_IntensityT2=I_const.*abs(signalInT2).^2;
% trapz(t,hcf_IntensityT2)
% trapz(f,hcf_IntensityF2)
% % scaling_Eloss=trapz(w,hcf_IntensityW)/trapz(w,hcf_IntensityW2)*0.4545;


%% Fourier Split Step Runge Kutta %% Fourier Split Step Runge Kutta 
n2_values=1e-24:0.5e-24:8e-24;
n2_length=length(n2_values);
Beta2_values=-3e-30:-0.5e-30:-6e-30;
Beta2_length=length(Beta2_values);

AbsDiff=zeros(Beta2_length,n2_length);
indexm=0;
indexl=0;
for m=n2_values 
    indexm=indexm+1; 
    for l=Beta2_values
        indexl=indexl+1;
        [prop_fieldT,prop_fieldF]=do_FourierSplitStep(f,signalInT,z,f0,l,m,I_const,alpha,'true');
        % prop_IntensityT=get_Intensity(prop_fieldT,I_const);
        prop_IntensityF=get_Intensity(prop_fieldF,I_const);
        
        
      
        AbsDiff(indexl,indexm)=calc_AbsoluteDifference(prop_IntensityF(408:618),hcf_IntensityF2(408:618));
    end
    indexl=0;
end 
% [prop_fieldT,prop_fieldF]=do_FourierSplitStep(f,signalInT,z,f0,Beta2,n2,I_const,alpha);
% prop_IntensityT=get_Intensity(prop_fieldT,I_const);
% prop_IntensityF=get_Intensity(prop_fieldF,I_const);

%% Convert back to Wavelength    
[lambda,prop_IntensityL]=convert_data(f,prop_IntensityF,'back');
[~,prop_IntensityL2]=convert_data(f,hcf_IntensityF2,'back');
    
%% Calculate absolute difference
absolute_diffenrece=calc_AbsoluteDifference(hcf_IntensityF2,prop_IntensityF);
%% Energy Conservation
% Input to Output
% E_cons=trapz(t,hcf_IntensityT)/trapz(t,hcf_IntensityT2);
% E_cons=trapz(f,hcf_IntensityF)/trapz(f,hcf_IntensityF2);

% Calculated to Input
% E_cons=trapz(f,hcf_IntensityF)/trapz(f,prop_IntensityF);

% Calculated to Output
% E_cons=trapz(t,hcf_IntensityT2)/trapz(t,prop_IntensityT);
% E_cons=trapz(f,hcf_IntensityF2)/trapz(f,prop_IntensityF);

%% Plot 

% plot(n2_values,AbsDiff)
figure; surf(n2_values,Beta2_values,AbsDiff./max(max(AbsDiff)));
ylabel('Beta2 [s^2/m]')
xlabel('n2 [m^2/W]')

% figure; imagesc(z,w,matrixw);
% title('In freq');
% % ylabel('w spectral beam spread [Thz]')
% % xlabel('z propagation direction [m]')
% 
% figure; imagesc(z,t,matrixt);
% title('In time');
% ylabel('t temporal beam spread [fs]')
% xlabel('z propagation direction [m]')
% % 
% figure;
figure;
subplot(2,1,1)
hold on
plot(lambda(400:1000).*1e9,[prop_IntensityL(400:1000);prop_IntensityL2(400:1000)])
legend('Theory','Data')
subplot(2,1,2)
hold on
plot(f,[prop_IntensityF;hcf_IntensityF2])
legend('Theory','Data','w0')

% figure; 
% plot(w,matrixw(:,end)./max(matrixw(:,end)));
% hold on
% plot(w,signalInW_2./max(signalInW_2));
% hold off
% figure;
% for s=1:1:length(z)
%     
%     
% 
%    plot(matrixw(:,s));
%    
%    axis([0 2000 0 max(max(matrixw))]);    
%     
%     pause(0.001);
% end

% plot(y_dataW)
% hold on
% plot(y_dataW2)
% hold off


