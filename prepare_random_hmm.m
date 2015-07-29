function [ dummy_hmm ] = prepare_random_hmm( F3, K )
    n_species = length( F3 );
    take_formants = 1;
    max_f = 0;
    min_f = 0;
    len = 0;
    %Calculate range
    for fmts = F3
        tmp_f = fmts{1}(:, [1:take_formants]);

        tmp_max = max( tmp_f );
        if max( tmp_f ) > max_f
            max_f = tmp_max;
        end

        tmp_min = min( tmp_f );
        if min( tmp_f ) < min_f 
            min_f = tmp_min;
        end

        len = len + length( tmp_f );
    end

    len = round( len / n_species );

    display('Training HMMs...');
    dummy_hmm.type = 'hmm';
    dummy_hmm.hmm = cell(1, n_species);
    dummy_dec = zeros(n_species, K);
    for i=1:n_species
        i
        hmm.K=K;
        data = get_random_formants( max_f, min_f, len );
        hmm = hmminit( data, hmm, 'full' );
        m = length( data );
        dummy_hmm.hmm{i} = hmmtrain( data, m, hmm );

        V = hmmdecode( data, length( data ), dummy_hmm.hmm{i} );
        for j=1:K
            dummy_dec(i, j) = sum(V.q_star == j);
        end

        [s ind] = sort(dummy_dec(i, :) , 'descend');
        dummy_hmm.hmm{i} = reorder_hmm_states( dummy_hmm.hmm{i}, ind );
        dummy_hmm.hmm{i}.data.occupancy = get_occupancy( dummy_hmm.hmm{i} );
        dummy_dec(i, :) = dummy_dec(i, ind);
        clear hmm;
    end

    
end