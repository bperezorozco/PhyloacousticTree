function [  ] = plot_subtree_genera( L, genus, id_folders )
%PLOT_SUBTREE_GENERA Summary of this function goes here
%   Detailed explanation goes here

index = find_ancestor_genus( genus, id_folders, L );
get_subtree(index, L, true, id_folders);
end

