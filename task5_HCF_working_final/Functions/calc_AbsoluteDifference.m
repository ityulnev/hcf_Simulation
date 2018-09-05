%Calculates the difference of absolutes between two datasets 
function [absolute_difference]=calc_AbsoluteDifference(data1,data2)
absolute_difference=sum(abs(data1-data2));
end