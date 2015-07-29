clusters = cell(2 * n_species - 1, 1);
clusters(1:n_species) =  id_folders(:) ;


for i=1:length(L)
    l = L(i, 1);
    r = L(i, 2);
    i
    clusters{i + n_species} = [ clusters(l); clusters(r) ];
    i + n_species
    clusters{i + n_species} 
end

figure;
plot( L(:, 3), (n_species-1):-1:1)


L2 = linkage(squareform(randn(1, n_species * (n_species - 1)/2)), 'average');
figure;
plot( L2(:, 3), (n_species-1):-1:1 )