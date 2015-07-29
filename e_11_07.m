display('Calculating distances between pairs of HMMs...');
d = zeros(n_species*(n_species-1)/2, K);
k = 1;
for i=1:n_species
    for j=i+1:n_species
        d(k, :) = distance_hmm( trained_hmm{i}, trained_hmm{j}, [1:K], 'partition' );
        k = k + 1;
    end
end

display('Performing hierarchical clustering...');
for take=1:10
    h = figure;
    L = linkage( squareform(d(:, take)), 'average' );
    dendrogram(L, n_species, 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3)));
    title(sprintf('Phylloacoustic tree using HMM+Dirichlet SKLD+take=%d and only birdsong', take));
end

display('Finished.');
clear i j f s perm hmm m V s ind dec L tmp