function [ X ] = unfold_cell_array( Cells )
%UNFOLD_CELL_ARRAY Summary of this function goes here
%   Detailed explanation goes here

cur = 1;
len = length(Cells);

for i=1:len
    M = Cells{i};
    [ r c ] = size( M );
    
    X(cur:(cur+r-1), :) = M;
    
    cur = cur + r;
end

end

