function [ sub order leaves ] = get_subtree( id_node, L, plot, labels )
%GET_SUBTREE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    plot = false;
end
if nargin == 3
    display('Required: @labels');
    return;
end
[r c] = size(L);
if id_node <= r+1
    %display('Singleton cluster');
    sub = [];
    order = [];
    leaves = 1;
    return;
end

order = [];
sub = span_subtree( id_node, L );
[sub order] = rename_subtree( sub, r + 1 );

leaves = length(order);
%display(sprintf('Number of leaves: %i', length(order)));
if plot
    figure;
    dendrogram(sub, length(order), 'Labels', labels(order), 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3))); 
    title(sprintf('Subtree for index: %i', id_node));
end

end

