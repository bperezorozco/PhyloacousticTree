K = 10;
n = length(trained_hmm);
d = zeros(n*(n-1)/2, K);
k = 1;
for i=1:n_species
    for j=i+1:n_species
        d(k, :) = distance_hmm(trained_hmm{i}, trained_hmm{j}, [1:K], 'partition');
        %d(k, :) = distance_hmm(trained_hmm{i}, trained_hmm{j}, ones(1, K), ones(1, K), [1:K], 'gauss');
        k = k + 1;
    end
end

for take=1:10
    h = figure;
    L = linkage( squareform(d(:, take)), 'average' );
    dendrogram(L, n_species, 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3)));
    title(sprintf('Phylloacoustic tree using HMM+Dirichlet SKLD+take=%d and only birdsong', take));
end
%title('Phylloacoustic tree using occupancy sorting+Dirichlet SKLD');

% A = zeros(length(L)*2+1, length(L)*2+1);
% count = length(L) + 2;
% for i=1:length(L)
%     A( L(i, 1),  count ) = 1;
%     A( L(i, 2),  count ) = 1;
%     count = count + 1;
% end
% A = A + A';
% 
% D = all_shortest_paths(sparse(A));
% t=5;
% k=1;
% for i=1:79
%     for j=1:79
%         if D(i,j) < t && D(i,j) > 0
%             s{k} = sprintf('(%s, %s)', folders5{i}, folders5{j});
%             display(s{k});
%             k=k+1;
%         end
%     end
% end