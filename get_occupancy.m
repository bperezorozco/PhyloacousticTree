function [ norm_dec ] = get_occupancy( trained_hmm )
%GET_OCCUPANCY Summary of this function goes here
%   Detailed explanation goes here
K = length(trained_hmm.state);
hmm_data = trained_hmm.data.Xtrain;

dec = zeros(1, K);
V = hmmdecode(hmm_data, length(hmm_data), trained_hmm);
for j=1:K
    dec(j) = sum(V.q_star == j);
end

norm_dec = dec / sum(dec);
end

