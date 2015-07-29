function [ V ] = train_hmm_BIC( seq, max_K )
%TRAIN_HMM_BIC Summary of this function goes here
%   Detailed explanation goes here


ll = zeros(1, max_K);
for i=1:max_K
    hmm = struct('K', i);
    hmm=hmminit( seq, hmm, 'full');
    [m , ~] = size( seq );
    hmm.train.rdisplay = false;
    
    hmm=hmmtrain( seq, m, hmm );
    
    [v, ~, l] = hmmdecode( seq, length( seq ), hmm );
    V(i) = v;
    nn(i) = numel( unique(v.q_star) );
    ll(i) = l;
end
nn
ll
end

