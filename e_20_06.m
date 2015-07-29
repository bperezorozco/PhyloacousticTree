clear d dec dec_dummy dec_dummy2 folder hmm hm_dummy i ind j m v V

folder = 5;
n_dummy = 5;
data = [];
K = 10;
for j=1:15
    data = [data; formants5(files5{(folder-1)*15+j})];
end

hmm_dummy = cell(1, n_dummy);

for i=1:n_dummy
    clear hmm;
    hmm.K=K;
    hmm=hmminit(data, hmm, 'full');
    m = length( data );
    hmm=hmmtrain(data, m, hmm);
    hmm_dummy{i} = hmm;
end

dec = zeros(n_dummy, K);
v = cell(1, n_dummy);

for i=1:n_dummy
    V = hmmdecode(data, length(data), hmm_dummy{i});
    
    for j=1:10
        dec(i, j) = sum(V.q_star == j);
    end
    
    [s ind] = sort(dec(i, :) , 'descend');
    %figure;
    %plot(V.q_star);
    v{i} = reorder_hmm_states(hmm_dummy{i}, ind, K);
    dec(i, :) = dec(i, ind);
end
%%
xmesh = 1:2^12-1;
Hhellinger = @(x,y)((1/sqrt(2))*sqrt(sum(bsxfun(@minus,sqrt(x),sqrt(y)).^2,2)));
Hskld = @(p,q)(sum(bsxfun(@times,bsxfun(@minus,p,q),log(bsxfun(@rdivide,p,q))),2));
pdf_gmm = zeros(length(trained_hmm), length(xmesh));
for i=1:length(trained_hmm)
    pdf_gmm(i, :) = gmmprob( trained_hmm{i}.mix, xmesh' )';
end

figure;
L = linkage( pdist(pdf_gmm, @(Xi, Xj)Hskld(Xi, Xj)), 'weighted' );
dendrogram(L, 79, 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3)));
title('Phylloacoustic tree using GMM+SKLD');
eval_phylo_tree(L, id_folders, meta)

figure;
L = linkage( pdist(pdf_gmm, @(Xi, Xj)Hhellinger(Xi, Xj)), 'weighted' );
dendrogram(L, 79, 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.25*max(L(:, 3)));
title('Phylloacoustic tree using GMM+Hellinger distance');
eval_phylo_tree(L, id_folders, meta)
% figure;
% hold on;
% plot(dec(1, :), 'b');
% plot(dec(2, :), 'r');
% plot(dec(3, :), 'g');
% plot(dec(4, :), 'k');
% plot(dec(5, :), 'c');