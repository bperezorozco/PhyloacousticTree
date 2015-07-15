n_files = 15;
K = 10;

hmm_data = cell(1, length(folders5));
for i=1:length(folders5)
x = [];
    for j=1:n_files
        x = [x; formants5(files5{(i-1)*n_files+j})];
    end
    hmm_data{i} = x;
end

trained_hmm = cell(1, length(hmm_data));
dec = zeros(length(hmm_data), K);
for i=1:length(hmm_data)
    clear hmm;
    hmm.K=K;
    hmm=hmminit(hmm_data{i}, hmm, 'full');
    m = length( hmm_data{i} );
    trained_hmm{i}=hmmtrain(hmm_data{i}, m, hmm);
    
    V = hmmdecode(hmm_data{i}, length(hmm_data{i}), trained_hmm{i});
    for j=1:K
        dec(i, j) = sum(V.q_star == j);
    end
    
    [s ind] = sort(dec(i, :) , 'descend');
    trained_hmm{i} = reorder_hmm_states(trained_hmm{i}, ind, K);
    dec(i, :) = dec(i, ind);
end