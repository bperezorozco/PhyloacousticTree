function [ s ] = get_random_TIER_file( filenames, folder_path )
%GET_RANDOM_TIER_FILE Summary of this function goes here
%   Detailed explanation goes here

s = strcat( folder_path, filenames{randi( length(filenames) )} );

end

