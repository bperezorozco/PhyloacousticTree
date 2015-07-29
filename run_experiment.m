function [ L res e ] = run_experiment( data, id_experiment, id_folders, params, meta, genera, dummy )
%RUN_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here
n_species = length( id_folders );
n_exp = 5;
if ~isfield(params, 'ltype')
    params.ltype = 'average';
end
if ~isfield(params, 'K')
    params.K = 10;
end
if ~isfield(params, 'xmesh')
    xmesh = [1:8192];
end

%get distance matrix

%Calculate distance between pdfs ( experiments between 0 and 9 )
if id_experiment < 10
    return;
%Calculate distance between GMMs ( experiments between 10 and 19 )
elseif id_experiment< 20 && strcmp(data.type, 'hmm') 
    xmesh = params.xmesh;
    pdf_gmm = zeros(length( data.hmm ), length(xmesh));
    pdf_gmm_random = zeros(length( data.hmm ), length(xmesh));
    trained_hmm = data.hmm;
    
    for i=1:length(trained_hmm)
        pdf_gmm(i, :) = gmmprob( trained_hmm{i}.mix, xmesh' )';
        pdf_gmm_random(i, :) = gmmprob( dummy.hmm{i}.mix, xmesh' )';
    end
    

    if mod( id_experiment, 10 ) == 1 && isfield( params, 'skld' ) && isfield( data, 'hmm' )
        L = linkage( pdist(pdf_gmm, @(Xi, Xj)params.skld(Xi, Xj)), params.ltype );
        L2 = linkage( pdist(pdf_gmm_random, @(Xi, Xj)params.skld(Xi, Xj)), params.ltype );
        titles = 'HMM GMM SKLD';
        titles_random = 'Random HMM GMM SKLD';
    elseif mod( id_experiment, 10) == 2 && isfield( params, 'hellinger' )
        pdist(pdf_gmm, @(Xi, Xj)params.hellinger(Xi, Xj))
        L = linkage( pdist(pdf_gmm, @(Xi, Xj)params.hellinger(Xi, Xj)), params.ltype )
        pdist(pdf_gmm_random, @(Xi, Xj)params.hellinger(Xi, Xj))
        L2 = linkage( pdist(pdf_gmm_random, @(Xi, Xj)params.hellinger(Xi, Xj)), params.ltype )
        titles = 'HMM GMM Hellinger';
        titles_random = 'Random HMM GMM Hellinger';
        figure;
        histogram(pdist(pdf_gmm, @(Xi, Xj)params.hellinger(Xi, Xj)))
        figure;
        histogram(pdist(pdf_gmm_random, @(Xi, Xj)params.hellinger(Xi, Xj)))
    else
        display('Experiment not supported yet. params.skld and params.hellinger must exist');
        return;
    end
    
%Calculate distance between HMMs ( experiments greater than 20 )
elseif id_experiment >= 20 && strcmp(data.type, 'hmm') && isfield( data, 'hmm' )
    if ~( isfield(params, 'K') )
        display('Missing param K in @params');
        return;
    end
    
    K = params.K;
    
    switch mod( id_experiment, 10 )
        case 1
            params.sort = 'mixture';
        case 2
            params.sort = 'weighted';
        case 3
            params.sort = 'occupancy';
        case 4
            params.sort = 'partition';
        otherwise
            params.sort = 'partition';
    end
    
    trained_hmm = data.hmm;
    params.take = [1:K];
    
    display('Calculating distances between pairs of HMMs...');
    d = zeros(n_species*(n_species-1)/2, params.K);
    random_d = zeros(n_species*(n_species-1)/2, params.K);
    k = 1;
    for i=1:n_species
        for j=i+1:n_species
            d(k, :) = distance_hmm( trained_hmm{i}, trained_hmm{j}, params );
            random_d(k, :) = distance_hmm( dummy.hmm{i}, dummy.hmm{j}, params );
            k = k + 1;
        end
    end


    display('Performing hierarchical clustering...');
    min_val = 10000;
    min_random = 10000;
    min_t = 1;
    min_t_random = 1;
    for take=1:K
        Ltmp = linkage( squareform(d(:, take)), 'average' );
        e = eval_phylo_tree2( Ltmp, genera, id_folders );
        
        if e < min_val
            L = Ltmp;
            min_val = e;
            min_t = take;
        end
        
        Ltmp = linkage( squareform(random_d(:, take)), 'average' );
        e = eval_phylo_tree2( Ltmp, genera, id_folders );
        if e < min_random
            L2 = Ltmp;
            min_random = e;
            min_t_random = take;
        end
    end
    
    titles = sprintf('Best result for HMM distance using %i states', min_t);
    titles_random = sprintf('Best result for HMM distance using %i states', min_t_random);
else
    display('Check DATA.TYPE is: gmm, hmm');
    return;
end

%show dendrogram
figure;
dendrogram(L, n_species, 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3)));
title(titles);

figure;
dendrogram(L2, n_species, 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3)));
title(titles_random);

max_distance = max( [ L(:, 3); L2(:, 3) ] );
n_points = 10000;

mesh = 0:(1/n_points):max_distance;
figure;
hold on;
plot(mesh, plot_threshold_clusters( L, max_distance, n_points ));
plot(mesh, plot_threshold_clusters( L2, max_distance, n_points ), 'r');
title('Number of clusters vs. lifetime');
xlabel('Lifetime');
ylabel('Clusters');

%run test of tree
res = analyse_linkage( L, genera, id_folders, n_exp );
display('Measure for measured tree');
e = eval_phylo_tree2( L, genera, id_folders )
display('Measure for random tree');
e2 = eval_phylo_tree2( L2, genera, id_folders )
end

