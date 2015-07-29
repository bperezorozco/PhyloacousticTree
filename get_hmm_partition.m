function [ p q minn ] = get_hmm_partition( D, symmetric )
%GET_HMM_PARTITION Summary of this function goes here
%   Detailed explanation goes here

K = length(D);

marked = zeros(K);
tmp_p = zeros(1, K);
tmp_q = zeros(K, 1);
p = [];
q = [];
D(D==0) = Inf;

dim = K;
for k=1:K
    tmp = reshape(D(~marked), dim, dim);
    minn(k) = min(tmp(:));
    [i j] = find( D == minn(k) );
    
    %allows to have repeated elements
    s = length(find(minn == minn(k)));
    
    i = i(s);
    j = j(s);
    
    %display(sprintf('%d with %d', i, j));
    
    
    tmp_p(j) = 1;
    tmp_q(i) = 1;
    p = [p i];
    q = [q j];
    
    marked = marked | repmat(tmp_p, K, 1) | repmat(tmp_q, 1, K);
    dim = dim - 1;
    
end

end

