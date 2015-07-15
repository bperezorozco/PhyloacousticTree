function [ hmm ] = reorder_hmm_states( hmm, indices )
%REORDER_HMM_STATES Summary of this function goes here
%   Detailed explanation goes here
K = length( hmm.state );

tmp = 0;
hmm.Dir_alpha = hmm.Dir_alpha(indices);

hmm.Dir2d_alpha = hmm.Dir2d_alpha(indices, :);
hmm.Dir2d_alpha = hmm.Dir2d_alpha(:, indices);

hmm.P = hmm.P(indices, :);
hmm.P = hmm.P(:, indices);

hmm.prior.Dir_alpha = hmm.prior.Dir_alpha(indices);

hmm.prior.Dir2d_alpha = hmm.prior.Dir2d_alpha(indices, :);
hmm.prior.Dir2d_alpha = hmm.prior.Dir2d_alpha(:, indices);

hmm.state = hmm.state(indices);

hmm.train.Gamma = hmm.train.Gamma(:, indices);
hmm.train.Gammasum = hmm.train.Gammasum(indices);
hmm.train.Xi = hmm.train.Xi(:, indices, :);
hmm.train.Xi = hmm.train.Xi(:, :, indices);

hmm.mix.priors = hmm.mix.priors(indices);
hmm.mix.centres = hmm.mix.centres(indices);
hmm.mix.covars = hmm.mix.covars(:, :, indices);

end

