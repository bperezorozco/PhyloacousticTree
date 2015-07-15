function [ d ] = skld( P, Q, X )
%SKLD Symmetric Kullback-Leibler divergence
%   Detailed explanation goes here

L = zeros(length(P), 1);
perform = bsxfun( @and, P, Q );
p = P(perform);
q = Q(perform);
d = sum( ( p - q ) .* log( p ./ q ), 2 );

end

