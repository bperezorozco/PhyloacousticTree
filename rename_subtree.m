function [ L order ] = rename_subtree( L, limit )
%RENAME_SUBTREE Summary of this function goes here
%   Detailed explanation goes here

[n c] = size(L);

query = (L(1:2*n) <= limit);
order = L(query);
actual_n = length(order);
[order indices] = sort(L(1:2*n), 'ascend');
L( indices ) = 1:(2*n);
order = order(1:actual_n);
end

