function [f] = mel2freq (m)
% f = mel2freq (m)
% compute frequency from mel value
f = 700*((10.^(m ./2595)) -1);
