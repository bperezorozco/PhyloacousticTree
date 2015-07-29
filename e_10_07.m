%% 

%load_birdsong;
load_metadata;
id_folders = unique(folders);
perm = randperm( length(files) );
j = 1;
r = 2;
K = 10;
max_secs = 20;
max_fs = 44100;
n_formants = 3;
n_species = length(id_folders);
timeframe = 40;
F3 = cell(1, n_species);
B3 = cell(1, n_species);
display('Preparing data and calculating formants...');
for i=1:length(perm)
    %s = strcat('../../dataset/', folders{perm(i)}, '/', files{perm(i)}, '.wav');
    i
    s = strcat('F:/dissertation_MAT/dataset/', folders{perm(i)}, '/', files{perm(i)}, '.wav');
    
    if ~(exist(s, 'file') == 2)
        display('File not found:');
        files{perm(i)}
        t{j} = files{perm(i)};
        j = j + 1;
        continue;
    end
    
    [x fs] = audioread( strcat('F:/dissertation_MAT/dataset/', folders{perm(i)}, '/', files{perm(i)}, '.wav') );
    
    if fs ~= max_fs
        display('File not sampled at 44.1 KHz');
        files{perm(i)}
        continue;
    end
    
    [tmp tmp2] = formants_from_mfcc( resample(x, 1, r), max_fs/r, n_formants, timeframe );
    idf = find(not(cellfun('isempty', strfind(id_folders, folders{ perm(i) }))));
    F3{idf} = [F3{idf}; tmp];
    %B3{idf} = [B3{idf}; tmp2];
    pause;
end

%%

display('Training HMMs...');
trained_hmm = cell(1, n_species);
dec = zeros(n_species, K);
for i=1:n_species
    i
    hmm.K=K;
    data = F3{i}(:, 1);
    hmm = hmminit( data, hmm, 'full' );
    m = length( data );
    trained_hmm{i} = hmmtrain( data, m, hmm );
    
    V = hmmdecode( data, length( data ), trained_hmm{i} );
    for j=1:K
        dec(i, j) = sum(V.q_star == j);
    end
    
    [s ind] = sort(dec(i, :) , 'descend');
    trained_hmm{i} = reorder_hmm_states( trained_hmm{i}, ind );
    trained_hmm{i}.data.occupancy = get_occupancy( trained_hmm{i} );
    dec(i, :) = dec(i, ind);
    clear hmm;
end
%% 
display('Calculating distances between pairs of HMMs...');
d = zeros(n_species*(n_species-1)/2, K);
k = 1;
for i=1:n_species
    for j=i+1:n_species
        d(k, :) = distance_hmm( trained_hmm{i}, trained_hmm{j}, [1:K], 'mixture' );
        k = k + 1;
    end
end


display('Performing hierarchical clustering...');
for take=1:10
    h = figure;
    L = linkage( squareform(d(:, take)), 'average' );
    eval_phylo_tree(L, id_folders, meta)
    dendrogram(L, n_species, 'Labels', id_folders, 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3)));
    title(sprintf('Phylloacoustic tree using HMM+Dirichlet SKLD+take=%d and only birdsong', take));
end

display('Finished.');
clear i j f s perm hmm m V s ind dec L tmp