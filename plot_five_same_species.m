while true
    n = randi(length(filenames));
    l = length( filenames{n} ) - 15;
    f1 = filenames{n};
    f2 = filenames{n+4};
    if f1(1:l) == f2(1:l)
        for i=1:5
            s{i} = strcat('../../', filenames{n + i - 1});
            %[F B] = formants_file( s{i}, 10, 5 );
            F = formants_from_mfcc_file( s{i}, 10, 5 );
            x{i} = F;
        end


        figure;
        title('Dist(F1, F2)');
        xlabel('F1 (Hz)')
        ylabel('F2 (Hz)');
        plot(x{1}(:, 1), x{1}(:, 2), 'b.', x{2}(:, 1), x{2}(:, 2), 'g.', x{3}(:, 1), x{3}(:, 2), 'r.', x{4}(:, 1), x{4}(:, 2), 'k.', x{5}(:, 1), x{5}(:, 2), 'm.');
        legend( s{1}, s{2}, s{3}, s{4}, s{5} );

        figure;
        hold on;
        title('F1 over time');
        plot(x{1}(:, 1), 'b.');
        plot(x{2}(:, 1), 'g.');
        plot(x{3}(:, 1), 'r.');
        plot(x{4}(:, 1), 'k.');
        plot(x{5}(:, 1), 'm.');
        legend( s{1}, s{2}, s{3}, s{4}, s{5} );
        
        figure;
        hold on;
        title('F2 over time');
        plot(x{1}(:, 2), 'b.');
        plot(x{2}(:, 2), 'g.');
        plot(x{3}(:, 2), 'r.');
        plot(x{4}(:, 2), 'k.');
        plot(x{5}(:, 2), 'm.');
        legend( s{1}, s{2}, s{3}, s{4}, s{5} );
        
        break;
    end
end
