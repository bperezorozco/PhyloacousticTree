function [ filtered new_bw ] = filter_formants( formants, threshold, take, bw )
%FILTER_FORMANTS Summary of this function goes here
%   Detailed explanation goes here

i = 1;
j = 1;
filtered = zeros(1, take);
new_bw = zeros(1, take);
while i <= length(formants) && j <= take
    if formants(i) > 0 
        %merge the formants separated by a margin of less than @threshold
        if i < length(formants) && abs(formants(i+1) - formants(i)) < threshold
            filtered(j) = filtered(j) + ( formants(i+1) + formants(i) ) / 2;
            new_bw(j) = new_bw(j) + ( bw(i+1) + bw(i) ) / 2;
            i = i + 1;
        else
            filtered(j) = filtered(j) + formants(i);
            new_bw(j) = new_bw(j) + bw(i);
        end
        j = j + 1;
    end
    
    i = i + 1;
end

%we don't want formants with this bw, so we just push them to be the last
%before sorting
new_bw( new_bw == 0 ) = 1000;

[new_bw ind] = sort(new_bw, 'ascend');
filtered = filtered(ind);
end

