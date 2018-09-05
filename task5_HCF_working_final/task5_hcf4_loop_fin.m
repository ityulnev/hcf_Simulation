close all
clear all
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
n0=calculate_refrindexNe(wavelength);
r_hcf=200e-6;%[m]
area_hcf=pi*(0.65*r_hcf)^2;%[m^2] Mode radius w~a*o.65
kfreq=n0*w0/const.c;%1/m
Fluence=Q_In/area_hcf;
alpha=0.7885;

n2=1.5e-24;%5.9/////3.6e-19*1e-4;%[m^2/W]  -> SPM
Beta2=-3.2215*(1e-30);%[s^2/m]           -> GVD

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

%% Calculate units based in center wavelength
weighted_middle=sum(x_data.*y_data)/sum(y_data);
wavelength=weighted_middle;
n0=calculate_refr_index(wavelength);
f0=const.c/wavelength;
w0=2*pi*f0;

%% Signal Amplitude
I_const=0.5*n0*const.c*const.eps0;

%% Prepare data for Fourier Split Step method
[x_dataF,y_dataF]=process_data(x_data,y_data,I_const);

%cut out data
cutrange=48:1048;%48:1048;%18:1100;%Middle for 800nm is 309:775  ****middle @Index 547-548 ***
f=x_dataF(cutrange);
y_dataFn=y_dataF(cutrange);
y_dataF=smooth_sides(y_dataFn,460,555);

df=f(2)-f(1);
t=linspace(-1/(2*df),1/(2*df),length(f));

scaling=sqrt(Fluence/(sum(get_Intensity(y_dataF,I_const))*df));
signalInF=y_dataF.*scaling;
signalInT=fftshift(ifft(ifftshift(signalInF))).*(length(f)*df);

hcf_IntensityT=I_const.*abs(signalInT).^2;
hcf_IntensityF=I_const.*abs(signalInF).^2;

%%Simulated gaussian input pulse
y_dataFtheo=exp(-(2*pi.*(f-f0)).^2.*tau^2./2);
scaling_theo=sqrt(Fluence/(sum(get_Intensity(y_dataFtheo,I_const))*df));
signalInFtheo=y_dataFtheo.*scaling_theo;
signalInTtheo=fftshift(ifft(ifftshift(signalInFtheo))).*(length(f)*df);
hcf_IntensityTtheo=I_const.*abs(signalInTtheo).^2;
hcf_IntensityFtheo=I_const.*abs(signalInFtheo).^2;

%% Real Data after HCF propagation
[x_dataF2,y_dataF2]=process_data(x_data2,y_data2,I_const);
y_dataF2=y_dataF2(cutrange);


scaling_Eloss=1/2.2;
scaling2=sqrt(scaling_Eloss*Fluence/(sum(get_Intensity(y_dataF2,I_const))*df));

signalInF2=y_dataF2.*scaling2;
hcf_IntensityF2=I_const.*abs(signalInF2).^2;
signalInT2=fftshift(ifft(ifftshift(signalInF2))).*(length(f)*df);%./sqrt(4*pi)
hcf_IntensityT2=I_const.*abs(signalInT2).^2;

%% Fourier Split Step Runge Kutta 
n2_values=2.2e-24:0.1e-24:2.6e-24;
n2_length=length(n2_values);
Beta2_values=-3e-30:-0.2e-30:-8e-30;
Beta2_length=length(Beta2_values);

AbsDiff=zeros(Beta2_length,n2_length);
indexm=0;
indexl=0;
for m=n2_values 
    indexm=indexm+1; 
    for l=Beta2_values
        indexl=indexl+1;
        [prop_fieldT,prop_fieldF]=do_FourierSplitStep(f,signalInTtheo,z,f0,l,m,I_const,alpha);
        % prop_IntensityT=get_Intensity(prop_fieldT,I_const);
        prop_IntensityF=get_Intensity(prop_fieldF,I_const);
        AbsDiff(indexl,indexm)=calc_AbsoluteDifference(prop_IntensityF(408:618),hcf_IntensityF2(408:618));
    end
    indexl=0;
end 
prop_IntensityT=get_Intensity(prop_fieldT,I_const);
prop_IntensityF=get_Intensity(prop_fieldF,I_const);

%% Convert back to Wavelength    
[lambda,hcf_IntensityL]=convert_data(f,prop_IntensityF,'back');
[~,hcf_IntensityL2]=convert_data(f,hcf_IntensityF2,'back');

%% Calculate absolute difference
absolute_diffenrece=calc_AbsoluteDifference(hcf_IntensityF2,prop_IntensityF);

%% Plot 

% plot(n2_values,AbsDiff)
figure; surf(n2_values,Beta2_values,AbsDiff./max(max(AbsDiff)));
ylabel('Beta2 [s^2/m]')
xlabel('n2 [m^2/W]')

% figure; imagesc(z,f,matrixf);
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
plot(lambda(400:1000).*1e9,[hcf_IntensityL(400:1000);hcf_IntensityL2(400:1000)])
ylabel('Intensity')
xlabel('Wavelength [nm]')
legend('Theory','Data')
subplot(2,1,2)
hold on
plot(f,[prop_IntensityF;hcf_IntensityF2])
ylabel('Intensity')
xlabel('Frequency [1/s]')
legend('Theory','Data','w0')

% figure; 
% plot(w,matrixw(:,end)./max(matrixw(:,end)));
% hold on
% plot(w,signalInW_2./max(signalInW_2));
% hold off
% figure;
% for s=1:1:length(z)
%    plot(matrixt(:,s));
%    
%    axis([0 2000 0 max(max(matrixt))]);    
%     
%     pause(0.001);
% end

% figure;
% energy_cons=zeros*(length(z));
% for s=1:1:length(z)   
%     energy_cons(1,s)=trapz(f,matrixf(:,s));       
% end
% test_conserved=energy_cons(1,1)/energy_cons(1,end);
% plot(z,energy_cons)
%%Tests
test_energy(t,f,prop_IntensityF,hcf_IntensityF,hcf_IntensityF2,prop_IntensityT,hcf_IntensityT,hcf_IntensityT2)

