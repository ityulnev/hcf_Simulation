%Converts Intensity in Units of Wavelength to Units of Energy or reverse
function [output_array_x,output_array_y]=convert_data(x_input_array,y_input_array,type)

%As there are no negative Intensities!
y_input_array(y_input_array<0)=0;

switch type
    case 'forward'
    output_array_x=(const.h*const.c)./x_input_array;
    output_array_y=(y_input_array.*(x_input_array.^2)./(const.h*const.c));
    case 'back'
    output_array_x=(const.c)./x_input_array;
    output_array_y=(y_input_array.*(const.h*const.c)./(output_array_x.^2));   
end


end