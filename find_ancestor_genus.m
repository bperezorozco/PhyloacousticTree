function [ index ] = find_ancestor_genus( genus, id_folders, L )
%FIND_ANCESTOR_GENUS Summary of this function goes here
%   Detailed explanation goes here

[ ind ] = find_indices_genus( genus, id_folders );
anc = cell(1, length(ind));
for i=1:length(ind)
    anc{i} = find_ancestors( ind(i), L );
end

index = find_first_common_ancestor(anc);
end

