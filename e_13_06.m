init;

s{1} = '../../dataset/Parus major/ParMaj00010.wav';
s{2} = '../../dataset/Parus major/ParMaj00024.wav';
s{3} = '../../dataset/Parus major/ParMaj00035.wav';
s{4} = '../../dataset/Parus major/ParMaj00052.wav';
s{5} = '../../dataset/Parus major/ParMaj00083.wav';
s{6} = '../../dataset/Parus major/ParMaj00140.wav';
s{7} = get_random_TIER_file( filenames, '../../' );
s{8} = get_random_TIER_file( filenames, '../../' );
s{9} = get_random_TIER_file( filenames, '../../' );
s{10} = get_random_TIER_file( filenames, '../../' );
s{11} = get_random_TIER_file( filenames, '../../' );
s{12} = get_random_TIER_file( filenames, '../../' );

estim_formants = true;
run = 1;

if run == 1
%%Experiment 1:
% K = 10
% Train on one example
% Use only the first formant
    if estim_formants
        clear F;
        for i=1:12
            tmp = formants_from_mfcc_file( s{i}, 6, 20 );
            F{i} = tmp(:, 1);
        end
    end
    
    clear hmm;
    hmm.K=10;
    hmm=hmminit(F{1}, hmm, 'full');
    [m n] = size( F{1} );
    hmm=hmmtrain(F{1}, m, hmm);

    for i=1:12
        [m n] = size( F{i} );
        [v, t, ll] = hmmdecode( F{i}, m, hmm);
        V{i} = v;
        LL(i) = ll;
    end

    [t, indices] = sort(LL);
    t
    s{indices}
elseif run == 2
%%Experiment 2:
% K = 4
% Each example is a frame in a signal
% Use first three formants

    if estim_formants
        clear F;
        for i=1:12
            tmp = formants_from_mfcc_file( s{i}, 6, 20 );
            F{i} = tmp(:, 1:3)';
        end
    end

    hmm.K=4;
    size( F{1} )
    hmm=hmminit(F{1}, hmm, 'full');
    hmm=hmmtrain( F{1}, 3, hmm);
    
    for i=1:12
        [m n] = size( F{i} );
        [v, t, ll] = hmmdecode( F{i}, m, hmm);
        V{i} = v;
        LL(i) = ll;
    end

    [t, indices] = sort(LL);
    s{indices}
end