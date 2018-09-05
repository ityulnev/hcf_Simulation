%Put in Electric Field in t and get back propagated field in w or t after
%distance z
%Implemented Effects: GVD, SPM and Energy Loss
function [prop_fieldT,prop_fieldF]=do_adaptiveFourierSplitStep(f,prop_fieldT,z,f0,Beta2,n2,I_const,alpha)

df=f(2)-f(1);
h=abs(z(2)-z(1));
w=f.*(2*pi);
w0=f0*2*pi;

% matrixt=[];
% matrixf=[];

prop_fieldT2=prop_fieldT;
error0=0.1;
d_prop=0;
    while (d_prop<=max(z))
                                      
            prop_fieldT =do_oneRKstep(w0,prop_fieldT,n2,I_const,h);                                                                                                                                                                                                                                                                                                                                         
            prop_fieldF=fftshift(fft(ifftshift(prop_fieldT)))./(length(f)*df);
            prop_fieldF=prop_fieldF.*exp((w-w0).^2.*Beta2/2*1i*h).*exp(-alpha/2*h); %GVD
            prop_fieldT=fftshift(ifft(ifftshift((prop_fieldF)))).*(length(f)*df);
    
            for m=1:1:2
            prop_fieldT2 =do_oneRKstep(w0,prop_fieldT2,n2,I_const,h/2);                                                                                                                                                                                                                                                                                                                                         
            prop_fieldF2=fftshift(fft(ifftshift(prop_fieldT2)))./(length(f)*df);
            prop_fieldF2=prop_fieldF2.*exp((w-w0).^2.*Beta2/2*1i*h/2).*exp(-alpha/2*h/2); %GVD
            prop_fieldT2=fftshift(ifft(ifftshift((prop_fieldF2)))).*(length(f)*df);
            end

            error=sum(abs(prop_fieldT-prop_fieldT2));
            if(error<error0)
                h=0.8*((error0/error)^(1/5));
            else
                h=0.8*((error0/error)^(1/4));
            end
            if(d_prop+h>max(z)&&h~=0)
                h=max(z)-d_prop;                             
            end
            d_prop=d_prop+h;  
%             matrixf=[matrixf;get_Intensity(prop_fieldF,I_const)];
%             matrixt=[matrixt;get_Intensity(prop_fieldT,I_const)];
    end

end

%% Runge Kutta outsourced functions

function [nextstep]=do_oneRKstep(w0,prop_fieldT,n2,I_const,h)
%calculate next step via Runge Kutta
       k1 = calcfunction(w0,prop_fieldT,n2,I_const);      
       k2 = calcfunction(w0,prop_fieldT+k1.*h./2,n2,I_const);        
       k3 = calcfunction(w0,prop_fieldT+k2.*h./2,n2,I_const);           
       k4 = calcfunction(w0,prop_fieldT+k3.*h,n2,I_const); 
       nextstep=prop_fieldT + h*(k1+2*k2+2*k3+k4)/6;
end

function [result]=calcfunction(w0,prop_fieldT,n2,Iconst)
%calculate diff. equation for Runge Kutta
result=1i*(n2*w0/const.c)*Iconst.*abs(prop_fieldT).^2.*prop_fieldT;
end