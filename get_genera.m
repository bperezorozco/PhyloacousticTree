function [ genera ] = get_genera( id_folders )
%GET_GENERA Summary of this function goes here
%   Detailed explanation goes here

genera = cell(1, length(id_folders));
for i=1:length(id_folders)
    s = strsplit( id_folders{i} );
    genera{i} = s{1};
end

genera = unique(genera);
end

