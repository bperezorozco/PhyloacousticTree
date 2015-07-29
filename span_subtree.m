function [ sub_L ] = span_subtree( id_node, L )
%SPAN_SUBTREE Summary of this function goes here
%   Detailed explanation goes here

[n c] = size(L);
sub_L = L( mod(id_node, n + 1), : );
done = false;

left=[];
right=[];
if ~(sub_L(2) <= n+1)
    left = span_subtree( sub_L(2), L );
end

if ~(sub_L(1) <= n+1)
    right = span_subtree( sub_L(1), L );
end
    
sub_L = [ left; right; sub_L ];
[s i] = sort( sub_L(:, 3) );

sub_L = sub_L(i, :);
end

