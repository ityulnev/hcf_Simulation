%Calculate linear refractive index in neon for wavelength 
%https://refractiveindex.info/?shelf=main&book=Ne&page=Bideau-Mehu
function [n0]=calculate_refrindexNe(wavelength)

%n0=1+0.029073880/(435.71376-wavelength^2); Cuthbertson's formula
n0=1+0.00128145/(184.661-(wavelength^(-2)))+0.0220486/(376.84-(wavelength^(-2)));

end