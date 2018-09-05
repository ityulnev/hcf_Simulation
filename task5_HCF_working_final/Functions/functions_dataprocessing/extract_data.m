%Takes a .fig File of a Plot and gets the Data of the X- and Y-Axis
function [x_In,y_In,x_Out,y_Out]=extract_data(figure_name,range,range2)

range=range(1,1):range(2,1);
range2=range2(1,1):range2(2,1);

% figure_name='HCFvsFund.fig';
open(figure_name);
figure_data=findobj(gcf,'type','line');

x_data=get(figure_data,'XData');
y_data=get(figure_data,'YData');
close(gcf);

x_In=x_data{2}(range);
y_In=y_data{2}(range);
x_Out=x_data{1}(range2);
y_Out=y_data{1}(range2);


end