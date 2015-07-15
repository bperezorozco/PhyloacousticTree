function [ E ] = count_links( S, k )
%COUNT_LINKS Summary of this function goes here
%   Detailed explanation goes here

n = length(S)
E = zeros( n );
for i=1:n
    r = S(:, i);
    [t indices] = sort(r(r>0));
    ix = indices(1:k);
    E(i, ix) = E(i, ix) + 1;
    E(ix, i) = E(ix, i) + 1;
end

end

