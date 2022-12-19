function output = replace_nan_with_avg(input,adjust) 
% function to replace NaNs in matrix with column averages. Used to give a
% peak latency for indivudal particiapnts in which not peak was found. 
%
% input  = N x M matrix of integers and Nans. N should be the number of
%          participants that the average will be computed from. M can be the
%          differnet conditions r regions 
% adjust = integer, used to adjust the peak latency to the correct value
%          based on the interval that was used to find the peak. Incorrect
%          because findpeaks returns index based on the interval we are given  
%
%
% output = N x M matrix of intergers same dimension as the input 

output            = input;
input_nan         = isnan(input);
avg_input         = round(nanmean(input));

for ii = 1:length(avg_input)
    temp = output(:,ii);
    temp(input_nan(:,ii)) = avg_input(ii);
    output(:,ii) = temp;
end
output = output + adjust; 

end