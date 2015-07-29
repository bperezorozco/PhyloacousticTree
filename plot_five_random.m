load filenames


for i=1:5
    s{i} = get_random_TIER_file(filenames, '../../');
    %[F B] = formants_file( s{i}, 10, 10 );
    F = formants_from_mfcc_file( s{i}, 10, 10 );
    x{i} = F;
    
    figure;
    hold on;
    plot(F(:, 1), 'b.')
    plot(F(:, 2), 'g.')
    plot(F(:, 3), 'r.')
    title( s{i} );
end


figure;
plot(x{1}(:, 1), x{1}(:, 2), 'b.', x{2}(:, 1), x{2}(:, 2), 'g.', x{3}(:, 1), x{3}(:, 2), 'r.', x{4}(:, 1), x{4}(:, 2), 'k.', x{5}(:, 1), x{5}(:, 2), 'm.');
legend( s{1}, s{2}, s{3}, s{4}, s{5} );

figure;
hold on;
plot(x{1}(:, 1), 'b.');
plot(x{2}(:, 1), 'g.');
plot(x{3}(:, 1), 'r.');
plot(x{4}(:, 1), 'k.');
plot(x{5}(:, 1), 'm.');
legend( s{1}, s{2}, s{3}, s{4}, s{5} );