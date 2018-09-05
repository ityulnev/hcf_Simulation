%Calculate linear refractive index in BK7 for wavelength 
%https://refractiveindex.info/?shelf=main&book=Ne&page=Bideau-Mehu
function [n0]=calculate_refrindexBK7(wavelength)

%n0=1+0.029073880/(435.71376-wavelength^2); Cuthbertson's formula
n0=sqrt(1+1.03961212*wavelength^2/(wavelength^2-0.00600069867)+0.231792344*wavelength^2/(wavelength^2-0.0200179144)+1.01046945*wavelength^2/(wavelength^2-103.560653));

end