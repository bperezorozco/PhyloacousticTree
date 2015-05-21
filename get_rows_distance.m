function [ D ] = get_rows_distance( X )
%GET_ROWS_DISTANCE Summary of this function goes here
%   Detailed explanation goes here

[m n] = size( X );
D = zeros(m, m);
temp = pdist(X);

k = 1;
for i=1:m
    D(i, i) = 0;
    for j=i+1:m
        D(i, j) = temp(k);
        D(j, i) = D(i, j);
        k = k + 1;
    end
end

end

