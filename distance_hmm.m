function [ d ] = distance_hmm( hmm_p, hmm_q, params )
%DISTANCE_HMM Summary of this function goes here
%   Detailed explanation goes here

take = params.take;
sort = params.sort;
K = length(hmm_p.state);

if nargin < 3
    sort = 'gauss';
    if nargin < 3
        take = K;
    end
end

if ~isfield( hmm_p.data, 'occupancy' )
    display('Warning: occupancy data not found in the first HMM. Calculating...');
    hmm_p.data.occupancy = get_occupancy( hmm_p );
    return;
end

if ~isfield( hmm_q.data, 'occupancy' )  
    display('Warning: occupancy data not found in the second HMM. Calculating...'); 
    hmm_q.data.occupancy = get_occupancy( hmm_q );
    return;
end

occ_p = hmm_p.data.occupancy;
occ_q = hmm_q.data.occupancy;

weights = zeros(1, K);
a = zeros(K);

if strcmp( sort, 'gauss' ) 
    for i=1:K
        minn = skld_unigauss(hmm_p.state(i).Mu, hmm_p.state(i).Cov, hmm_q.state(K).Mu, hmm_q.state(K).Cov);
        
        swap_ind = K;
        for j=i+1:K-1
            tmp = skld_unigauss(hmm_p.state(i).Mu, hmm_p.state(i).Cov, hmm_q.state(j).Mu, hmm_q.state(j).Cov);
            if tmp < minn
                minn = tmp;
                swap_ind = j;
            end
        end
        
        if minn == 0
            weights(i) = 0;
        else
            weights(i) = 1 / minn;
        end
        
        hmm_q = swap_hmm_states( hmm_q, i, swap_ind );
        occ_q([i swap_ind]) = occ_q([swap_ind i]);
    end
elseif strcmp(sort, 'partition') || strcmp(sort, 'mixture') || strcmp(sort, 'weighted')
    for i=1:K
        for j=1:K
            a(i, j) = skld_unigauss(hmm_p.state(i).Mu, hmm_p.state(i).Cov, hmm_q.state(j).Mu, hmm_q.state(j).Cov);
        end
    end
    [indp indq weights] = get_hmm_partition(a, false);
   
    weights( weights == 0 ) = 1;
    hmm_p = reorder_hmm_states( hmm_p, indp );
    hmm_q = reorder_hmm_states( hmm_q, indq );
    
    occ_p = occ_p(indp);
    occ_q = occ_q(indq);
end

%weights
%calculate dirichlet sklds ()
tmp = zeros(1, K);
if strcmp(sort, 'mixture')
    tmp = weights;
elseif strcmp(sort, 'weighted')
    d_i = occ_p .* prod(a);
    d_j = occ_q .* prod(a');
    tmp = d_i + d_j;
else
    for i=1:K
        %tmp(i) = skld_dirichlet( occ_p .* hmm_p.Dir2d_alpha(i, :), occ_q .* hmm_q.Dir2d_alpha(i, :) );
        tmp(i) = skld_dirichlet( hmm_p.Dir2d_alpha(i, :), hmm_q.Dir2d_alpha(i, :) );
    end
    
    tmp =  (occ_p + occ_q)/2 .* tmp;
end
%


% weights = 1 ./ weights;
% nc = sum(weights);
% if nc == 0
%     nc = 1;
% end
% tmp = (weights/sum(weights)) .* (tmp/norm(tmp));

d = zeros(1, length(take));
for t=1:length(take)
    %d(t) = sum(tmp);
    d(t) = sum( tmp(1:take(t)) );
end

end

