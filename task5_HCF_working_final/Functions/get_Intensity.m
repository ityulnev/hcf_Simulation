%Calculates Intensity from field 
function [intensity]=get_Intensity(field,I_const)
intensity=I_const.*abs(field).^2;
end