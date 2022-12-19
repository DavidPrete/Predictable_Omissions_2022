
function output = extract_means(input, peak_loc)
% input    = erp data in fieldtrip format after ft_timelockanalysis 
% peak_loc = 1x3 vector containing the ind of the peak latency for the
%            left, midline and right electrode regions 
%
%output    = 1x3 vector of the mean amplitude average around the peak
%            latencies provided by peak_loc. 


    %+/- 3 ensure the window of average is 20ms around the peak latency
    win_left  = peak_loc-2;
    win_right = peak_loc+2;
    
    output    = mean(input(1,win_left:win_right));

end