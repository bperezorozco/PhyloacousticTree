function [ hmm ] = swap_hmm_states( hmm, i, j )
%SWAP_HMM_STATES Summary of this function goes here
%   Detailed explanation goes here

tmp = 0;
hmm.Dir_alpha([i j]) = hmm.Dir_alpha([j i]);

hmm.Dir2d_alpha([i j], :) = hmm.Dir2d_alpha([j i], :);
hmm.Dir2d_alpha(:, [i j]) = hmm.Dir2d_alpha(:, [j i]);

hmm.P([i j], :) = hmm.P([j i], :);
hmm.P(:, [i j]) = hmm.P(:, [j i]);

hmm.prior.Dir_alpha([i j]) = hmm.prior.Dir_alpha([j i]);

hmm.prior.Dir2d_alpha([i j], :) = hmm.prior.Dir2d_alpha([j i], :);
hmm.prior.Dir2d_alpha(:, [i j]) = hmm.prior.Dir2d_alpha(:, [j i]);

hmm.state([i j]) = hmm.state([j i]);

hmm.train.Gamma(:, [i j]) = hmm.train.Gamma(:, [j i]);
hmm.train.Gammasum([i j]) = hmm.train.Gammasum([j i]);
hmm.train.Xi(:, [ i j ], :) = hmm.train.Xi(:, [ j i ], :);
hmm.train.Xi(:, :, [ i j ]) = hmm.train.Xi(:, :, [ j i ]);

hmm.mix.priors([i j]) = hmm.mix.priors([j i]);
hmm.mix.centres([i j]) = hmm.mix.centres([j i]);
hmm.mix.covars(:, :, [i j]) = hmm.mix.covars(:, :, [j i]);

% 
% new_hmm = hmm;
% tmp = hmm.Dir_alpha(i);
% new_hmm.Dir_alpha(i) = hmm.Dir_alpha( j ); new_hmm.Dir_alpha(j) = tmp;
% 
% tmp = hmm.Dir2d_alpha(i, :);
% new_hmm.Dir2d_alpha(i, :) = hmm.Dir2d_alpha(j, :); new_hmm.Dir2d_alpha(j, :) = tmp;
% tmp = new_hmm.Dir2d_alpha(:, i);
% new_hmm.Dir2d_alpha(:, i) = new_hmm.Dir2d_alpha(:, j); new_hmm.Dir2d_alpha(:, j) = tmp;
% 
% tmp = hmm.P(i, :);
% new_hmm.P(i, :) = hmm.P(j, :); new_hmm.P(j, :) = tmp;
% tmp = new_hmm.P(:, i);
% new_hmm.P(:, i) = new_hmm.P(:, j); new_hmm.P(:, j) = tmp;
% 
% tmp = 0;
% new_hmm = hmm;
% tmp = hmm.prior.Dir_alpha(i);
% new_hmm.prior.Dir_alpha(i) = hmm.prior.Dir_alpha( j ); new_hmm.prior.Dir_alpha(j) = tmp;
% 
% tmp = hmm.prior.Dir2d_alpha(i, :);
% new_hmm.prior.Dir2d_alpha(i, :) = hmm.prior.Dir2d_alpha(j, :); new_hmm.prior.Dir2d_alpha(j, :) = tmp;
% tmp = new_hmm.prior.Dir2d_alpha(:, i);
% new_hmm.prior.Dir2d_alpha(:, i) = new_hmm.prior.Dir2d_alpha(:, j); new_hmm.prior.Dir2d_alpha(:, j) = tmp;
% 
% tmp = hmm.state( i );
% new_hmm.state(i) = hmm.state(j); new_hmm.state(j) = tmp;
% 
% tmp = hmm.mix.priors( i );
% new_hmm.mix.priors(i) = hmm.mix.priors(j); new_hmm.mix.priors(j) = tmp;
% 
% tmp = hmm.mix.centres( i );
% new_hmm.mix.centres(i) = hmm.mix.centres(j); new_hmm.mix.centres(j) = tmp;
% 
% tmp = hmm.mix.covars(:, :, i);
% new_hmm.mix.covars(:, :, i) = hmm.mix.covars(:, :, j); new_hmm.mix.covars(:, :, j) = tmp;
% 
% tmp = hmm.train.Gamma(:, i);
% new_hmm.train.Gamma(:, i) = hmm.train.Gamma(:, j); new_hmm.train.Gamma(:, j) = tmp;
% 
% tmp = hmm.train.Gammasum(:, i);
% new_hmm.train.Gammasum(:, i) = hmm.train.Gammasum(:, j); new_hmm.train.Gammasum(:, j) = tmp;
% 
% tmp = hmm.train.Xi(:, :, i);
% new_hmm.train.Xi(:, :, i) = hmm.train.Xi(:, :, j); new_hmm.train.Xi(:, :, j) = tmp;
% 
% tmp = new_hmm.train.Xi(:, i, :);
% new_hmm.train.Xi(:, i, :) = new_hmm.train.Xi(:, j, :); new_hmm.train.Xi(:, j, :) = tmp;

end

