function [ d ] = skld_unigauss( m1, s1, m2, s2 )
%SKLD_UNIGAUSS Summary of this function goes here
%   Detailed explanation goes here

%l1 = 1 / s1;
%l2 = 1 / s2;

%d = ( kld_unigauss( m1, l1, m2, l2 ) + kld_unigauss( m2, l2, m1, l1 ) ) / 2;
d = ( gauss_kl( m1, m2, s1, s2 ) + gauss_kl( m2, m1, s2, s1 ) ) / 2;
end

function [ d ] = kld_unigauss( m1, l1, m2, l2 )

d = ( log( l2 / l1 ) - 1 + l1 / l2 + l1 * ( m2 - m1 ) ^ 2 ) / 2;

end
