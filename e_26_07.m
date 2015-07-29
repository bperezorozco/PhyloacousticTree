pdf_gmm = zeros(length(trained_hmm), length(xmesh));
for i=1:length(trained_hmm)
    pdf_gmm(i, :) = gmmprob( trained_hmm{i}.mix, xmesh' )';
end

figure;
L = linkage( pdist(pdf_gmm, @(Xi, Xj)Hskld(Xi, Xj)), 'weighted' );
dendrogram(L, length(id_folders), 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3)));
title('Phylloacoustic tree using KDE+SKLD');

figure;
L = linkage( pdist(pdf_gmm, @(Xi, Xj)Hhellinger(Xi, Xj)), 'weighted' );
dendrogram(L, length(id_folders), 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.45*max(L(:, 3)));
title('Phylloacoustic tree using KDE+Hellinger distance');

figure;
plot( L(:, 3), (n_species-1):-1:1)

L2 = linkage(squareform(randn(1, n_species * (n_species - 1)/2)), 'average');
figure;
plot( L2(:, 3), (n_species-1):-1:1 )
dendrogram(L2, length(id_folders), 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.45*max(L(:, 3)));