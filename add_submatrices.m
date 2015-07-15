function [ Conv ] = add_submatrices( M, step )
%CONVOLVE Summary of this function goes here
%   Detailed explanation goes here

[rM cM] = size( M );
rConv = rM / step + rem( rM, step );
cConv = cM / step + rem( cM, step );

tmp = zeros( rM + rem(rM, step), cM + rem(cM, step) );
tmp(1:rM, 1:cM) = tmp(1:rM, 1:cM) + M;

M = tmp;
Conv = zeros( rConv, cConv );

for i=0:rConv-1
    for j=0:cConv-1
        rbegin = i*step+1;
        rend = (i+1)*step;
        cbegin = j*step+1;
        cend = (j+1)*step;
        tmp = M(rbegin:rend, cbegin:cend);
        Conv(i+1,j+1) = sum( tmp(:) ); 
    end
end
