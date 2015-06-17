function [ H ] = hellinger( P, Q, X )
%HELLINGER Calculates the Hellinger distance between two pdfs
%   Detailed explanation goes here
    H = ( 1 / sqrt(2) ) * norm( sqrt(P) - sqrt(Q) );
end

