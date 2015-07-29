load filenames
for i=1:5
    s{i} = filenames{randi([1 length(filenames)])};
    x{i} = formants_file( strcat('../../', s{i}) );
end

plot3(x{1}(:, 1), x{1}(:, 2), x{1}(:, 3), 'b', x{2}(:, 1), x{2}(:, 2), x{2}(:, 3), 'g', x{3}(:, 1), x{3}(:, 2), x{3}(:, 3), 'r', x{4}(:, 1), x{4}(:, 2), x{4}(:, 3), 'k', x{5}(:, 1), x{5}(:, 2), x{5}(:, 3), 'm');
legend( s{1}, s{2}, s{3}, s{4}, s{5} );