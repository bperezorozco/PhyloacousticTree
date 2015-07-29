function [ e ] = eval_phylo_tree2( L, genera, id_folders )
%EVAL_PHYLO_TREE2 Summary of this function goes here
%   Detailed explanation goes here

e = 0;
for genusc=genera
    genus = genusc{1};
    card_g = length( find_indices_genus( genus, id_folders ) );
    index = find_ancestor_genus( genus, id_folders, L );
    [ sub order c_g ] = get_subtree( index, L, false, id_folders );
    e = e + log( c_g / card_g );
end

end

