n_exp = 5;

truth = zeros( 1, length(genera) );
expe = zeros( 1, length(genera) );
r = zeros( length(genera), n_exp );

L2 = linkage(squareform(randn(1, n_species * (n_species - 1)/2)), 'average');
for i=1:length(genera)
    genera{i}
    index =  find_ancestor_genus(genera{i}, id_folders, L);
    if index <= n_species
        continue;
    end
    
    [ tmp tmp2 leaves ] = get_subtree( index, L );
    truth(i) = leaves;
    
    [ tmp tmp2 leaves ] = get_subtree( index, L2 );
    expe(i) = leaves;
end

[truth; expe]'
