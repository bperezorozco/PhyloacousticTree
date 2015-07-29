function [ h ] = eval_phylo_tree( L, labels, meta )
%EVAL_PHYLO_TREE Summary of this function goes here
%   Detailed explanation goes here

%meta must have 4 arrays (in this order): species in english, species in
%latin, family and order

if ~isfield(meta, 'latin') || ~isfield(meta, 'fam') || ~isfield(meta, 'ord')
    display('Error: incorrect metadata supplied.');
    return;
end

%obtain genera
meta.genus = cell(1, length(meta.latin));
for i=1:length(meta.latin)
    sci = strsplit( meta.latin{i} );
    meta.genus{i} = sci{1};
end

n_labels = length( labels );
n_clusters = 2*n_labels - 1;
[rows cols] = size( L );
clusters =  repmat( struct('init',true) , n_clusters, 1 );
%each cluster has a subset of the supplied metadata.

%init the metadata for the individual clusters
for i=1:n_labels
    idf = find(not(cellfun('isempty', strfind(meta.latin, labels{i}))));
    clusters(i).genera = labels{i};
    clusters(i).fam = meta.fam{idf};
    clusters(i).ord = meta.ord{idf};
    clusters(i).val = 1;
    clusters(i).numel = 1;
    %display(strcat(labels{i}, ' belongs to  ', meta.fam{idf}, ' and  ', meta.ord{idf}));
end

for i=1:rows
    %merge clusters L(i, 1) and L(i, 2)
    %this means: union of their families, union of their orders, calculate
    %new index
    ind = n_labels + i;
    clusters(ind).genera = union( clusters(L(i, 1)).genera, clusters(L(i, 2)).genera);
    clusters(ind).fam = union( clusters(L(i, 1)).fam, clusters(L(i, 2)).fam);
    clusters(ind).ord = union( clusters(L(i, 1)).ord, clusters(L(i, 2)).ord);
    clusters(ind).numel = clusters(L(i, 1)).numel + clusters(L(i, 2)).numel;
    
    k =  clusters(ind).numel * length(clusters(ind).genera) * length(clusters(ind).fam) * length(clusters(ind).ord);
    weighted_sum = clusters(L(i, 1)).numel*clusters(L(i, 1)).val + clusters(L(i, 2)).numel*clusters(L(i, 2)).val;
    clusters(n_labels + i).val = ( weighted_sum ) / k;
end

h = -log( clusters(n_clusters).val );
end

