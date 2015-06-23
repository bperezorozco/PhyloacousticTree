function [ new_hmm ] = reorder_hmm_states( hmm, indices, K )
%REORDER_HMM_STATES Summary of this function goes here
%   Detailed explanation goes here

new_hmm = hmm;
new_hmm.Dir_alpha = hmm.Dir_alpha( indices );

new_hmm.Dir2d_alpha = hmm.Dir2d_alpha(indices, :);
new_hmm.Dir2d_alpha = new_hmm.Dir2d_alpha(:, indices);

new_hmm.P = hmm.P(indices, :);
new_hmm.P = new_hmm.P(:, indices);

new_hmm.state = hmm.state( indices );
new_hmm.mix.priors = hmm.mix.priors( indices );
new_hmm.mix.centres = hmm.mix.centres( indices );

new_hmm.mix.covars = hmm.mix.covars(:, :, indices);
new_hmm.train.Gamma = hmm.train.Gamma(:, indices);
new_hmm.train.Gammasum = hmm.train.Gammasum(:, indices);

new_hmm.train.Xi = hmm.train.Xi(:, :, indices);
new_hmm.train.Xi = new_hmm.train.Xi(:, indices, :);

end

