%Calculate the mode radius inside a step-index single mode fiber
%Marcuse Formula https://www.rp-photonics.com/mode_radius.html 
function [mode_radius]=calc_moderadius(wavelength,hcf_radius,n0,n0_fiber)

V=2*pi*hcf_radius*sqrt(n0^2-n0_fiber^2)/wavelength;
mode_radius=hcf_radius*(0.65+1.619/(V^(3/2))+2.879/(V^6));

end