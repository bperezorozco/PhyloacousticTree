function [ C ] = list_clusters( L, tags, n )
%LIST_CLUSTERS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    n = 10;
end

T = cluster(L, 'maxclust', n);
C = cell(1, n);
for i=1:n
    s = sprintf('Cluster %d', i);
    display(s);
    tags{T==i}
end

end

