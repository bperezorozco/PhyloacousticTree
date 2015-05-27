function [ F ] = folder_formants( folder )
%FOLDER_FORMANTS Summary of this function goes here
%   Detailed explanation goes here

i = 1;
n_bins = 10;
folder
s = strcat(folder, '/*.wav');
for file=dir( s )'
    file.name
    [x fs] = audioread( strcat(folder, '/', file.name) );
    F{i} = formants(x, fs);
    %F(i, :) = mean(formants(x, fs), 1);
    %t = formants(x, fs);
    %h = histogram(t(:), n_bins);
    %F(i, :) = h.Values / sum(h.Values);
    i = i + 1;
end

end

