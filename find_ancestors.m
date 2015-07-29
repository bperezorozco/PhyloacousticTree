function [ ancestors ] = find_ancestors( id, L )
%FIND_ANCESTORS Summary of this function goes here
%   Detailed explanation goes here

[n c] = size( L );
ancestors = [ id ];
row = mod( find(L == id), n );
while row ~= 0
    ancestor = (row+n+1);
    ancestors = [ ancestors ancestor ];
    row = mod( find(L == ancestor), n );
end

ancestors = [ ancestors 2*n+1 ];
end

%%
%obten numero de fila
% si el numero de fila 