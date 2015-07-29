function [ f ] = get_random_formants( max_f, min_f, k )
%GET_RANDOM_FORMANTS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
     k = 10000;
end

f = min_f + ( max_f - min_f ) .* rand( k, 1 );
end

