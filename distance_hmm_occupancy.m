function [ d ] = distance_hmm_occupancy( hmm_p, hmm_q )
%DISTANCE_HMM_OCCUPANCY Summary of this function goes here
%   Detailed explanation goes here
K = length( hmm_p.state );
d = 0;
for i=1:K
    d = d + skld_dirichlet( hmm_p.Dir2d_alpha(i, :), hmm_q.Dir2d_alpha(i, :) );
end

end

