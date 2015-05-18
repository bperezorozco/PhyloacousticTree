function [ F ] = folder_formants( folder )
%FOLDER_FORMANTS Summary of this function goes here
%   Detailed explanation goes here

i = 1;
for file=dir( strcat(folder, '/*.wav') )'
    [x fs] = audioread( strcat(folder, '/', file.name) );
    F{i} = formants(x, fs);
    mean(F{i}, 1)
    i = i + 1;
end

end

