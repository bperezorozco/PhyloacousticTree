function [ res ] = analyse_linkage( L, genera, id_folders, n_exp )
%ANALYSE_LINKAGE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    n_exp = 5;
end

n_species = length( id_folders )
truth = zeros( 1, length(genera) );
r = zeros( length(genera), n_exp );

for j=1:n_exp
    L2 = linkage(squareform(randn(1, n_species * (n_species - 1)/2)), 'average');
    
    
    for i=1:length(genera)
        index =  find_ancestor_genus(genera{i}, id_folders, L);
        index2 =  find_ancestor_genus(genera{i}, id_folders, L2);
        if index <= n_species
            truth(i) = 1;
        else
            [ tmp tmp2 leaves ] = get_subtree( index, L, false, id_folders );
            truth(i) = leaves;
        end
        
        if index2 <= n_species
            r(i, j) = 1;
        else
            [ tmp tmp2 leaves ] = get_subtree( index2, L2, false, id_folders );
            r(i, j) = leaves;
        end
        
    end
    
    
end
res = [truth' mean(r, 2)];
end

