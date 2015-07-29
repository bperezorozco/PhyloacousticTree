function [ ind ] = find_indices_genus( genus, id_folders )
%FIND_INDICES_GENUS Summary of this function goes here
%   Detailed explanation goes here

ind = find(not(cellfun('isempty', strfind(id_folders, genus))));

end

