%Put in Electric Field in t and get back propagated field in w or t after
%distance z
%Implemented Effects: GVD, SPM and Energy Loss via alpha
function [prop_fieldT,prop_fieldF,matrixt,matrixf]=do_FourierSplitStep(f,prop_fieldT,z,f0,Beta2,n2,I_const,alpha)

df=f(2)-f(1);
dz=abs(z(2)-z(1));
h=dz;
w=f.*(2*pi);
w0=f0*2*pi;

matrixt=zeros(length(f),length(z));
matrixf=zeros(length(f),length(z));

    for m=1:(length(z))
                 
            prop_fieldT=do_oneRKstep(w0,prop_fieldT,n2,I_const,h,w);                                                                                                                                                                                                                                                                                                                                         
            prop_fieldF=fftshift(fft(ifftshift(prop_fieldT)))./(length(f)*df);
            prop_fieldF=prop_fieldF.*exp((w-w0).^2.*Beta2/2*1i*dz).*exp(-alpha/2*dz); %GVD
            prop_fieldT=fftshift(ifft(ifftshift((prop_fieldF)))).*(length(f)*df);
            
            matrixf(:,m)=I_const.*abs(prop_fieldF).^2;
            matrixt(:,m)=I_const.*abs(prop_fieldT).^2;
    end
end

%% Runge Kutta outsourced functions

function [nextstep]=do_oneRKstep(w0,prop_fieldT,n2,I_const,h,w)
%calculate next step via Runge Kutta
       k1 = calcfunction(w0,prop_fieldT,n2,I_const,w);      
       k2 = calcfunction(w0,prop_fieldT+k1.*h./2,n2,I_const,w);        
       k3 = calcfunction(w0,prop_fieldT+k2.*h./2,n2,I_const,w);           
       k4 = calcfunction(w0,prop_fieldT+k3.*h,n2,I_const,w); 
       nextstep=prop_fieldT + h*(k1+2*k2+2*k3+k4)/6;
end

function [result]=calcfunction(w0,prop_fieldT,n2,I_const,w)
%calculate diff. equation for Runge Kutta
df=abs(w(2)-w(1))/(2*pi);
result=1i*(n2/const.c)*I_const.*abs(prop_fieldT).^2.*prop_fieldT;
result=w.*fftshift(fft(ifftshift(result)))./(length(w)*df);
result=fftshift(ifft(ifftshift((result)))).*(length(w)*df);
end