%Process data
function [x_data4,y_data4]=process_data(x_data,y_data,I_const)

%Convert from wavelength to energy
[x_data2,y_data2]=convert_data(x_data,y_data,'forward');
y_data2c=get_Efield(y_data2,I_const);
%Interpolate for evenly spaced Signal w,E(w)
[x_data3,y_data3]=even_data(x_data2,y_data2c);
%Unit coversion of x-axis  [J->1/s]
x_data3c=x_data3./const.h;
%Extending the W domain so that Simulation doesn't hit boundaries
[x_data4,y_data4]=extendToZero(x_data3c,y_data3);

%%
test_DataProcessing(x_data,y_data,x_data2,y_data2,y_data2c,x_data3,y_data3,x_data3c,x_data4,y_data4)

end