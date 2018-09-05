%Takes noisy data and smoothens left and right side of signal by setting
%noise to zero
function [y_data]=smooth_sides(y_data,signal_start,signal_end)
y_data(1:signal_start)=0;
y_data(signal_end:length(y_data))=0;
end

