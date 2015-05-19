function [ F ] = folder_formants( folder )
%FOLDER_FORMANTS Summary of this function goes here
%   Detailed explanation goes here

i = 1;
folder
s = strcat(folder, '/*.wav');
for file=dir( s )'
    file.name
    [x fs] = audioread( strcat(folder, '/', file.name) );
    F(i, :) = mean(formants(x, fs), 1);
    i = i + 1;
end

end

