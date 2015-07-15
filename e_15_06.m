init;
j=1;
k=1;
nfiles = 15;
min_length = 20;
display('Choosing files...');
for i=1:length(folders)
    s = strcat('../../dataset/', folders{i}, '/*.wav');
    d = dir( s )';
    if length(d) > min_length
        folders5{j} = folders{i};
        j = j+1;
        seq = randperm(length(d), nfiles);
        for h=1:nfiles
            files5{k} = strcat('dataset/', folders{i}, '/', d(seq(h)).name);
            k = k+1;
        end
    end
end

display('Obtaining formants...');
F1_5 = cell( 1, length( files5 ) );
for i=1:length(files5)
    F1_5{i} = formants(files5{i});
end
formants5 = containers.Map(files5, F1_5);

display('Estimating pdfs...');
pdf5 = cell( 1, length( files5 ) );
for j=1:length(files5)
    %[bandwidth,density,xmesh] = kde( formants( files5{j} ), 0, 8000 );
    [density xi bw] = ksdensity( formants( files5{j} ), [1:8192]  );
    pdf5{j}.bw = bw;
    pdf5{j}.density = density;
end
xmesh = xi;

display('Estimating average pdf per species...');
pdf_sp = zeros( length( folders5 ), length( pdf5{1}.density ) );
for j=1:length(folders5)
accum = 0;
    for k=1:nfiles
        accum = accum + pdf5{j*(nfiles-1)+k}.density;
    end
pdf_sp(j, :) = mean(accum / nfiles);
end

% display('Calculating SKLD and Hellinger distance matrices...');
Hhellinger = @(x,y)((1/sqrt(2))*sqrt(sum(bsxfun(@minus,sqrt(x),sqrt(y)).^2,2)));
Hskld = @(p,q)(sum(bsxfun(@times,bsxfun(@minus,p,q),log(bsxfun(@rdivide,p,q))),2));

% S = get_rows_distance(pdf5, Hskld);
% H = get_rows_distance(pdf5, Hhellinger);
% 
% 
% display('Calculating SKLD and Hellinger distance matrices per species...');
% S2 = get_rows_distance(pdf_sp, Hskld);
% H2 = get_rows_distance(pdf_sp, Hhellinger);

% display('Plotting K species pdfs at random...');
% K = 6;
% seq = randperm(length(folders5), K); 
% plot( xmesh, pdf_sp( seq, : ) ); 
% legend(folders5{seq});
% 
% figure;
% pcolor(exp(-H2));shading('flat')
% figure;
% pcolor(exp(-S2));shading('flat')
% clear bandwidth d density F1 filenames folders formants h i j k limit min_length nfiles s seq;
% display('Finished.');

figure;
L = linkage( pdist(pdf_sp, @(Xi, Xj)Hskld(Xi, Xj)), 'weighted' );
dendrogram(L, 79, 'Labels', folders5, 'Orientation', 'right', 'ColorThreshold', 0.3*max(L(:, 3)));
title('Phylloacoustic tree using KDE+SKLD');

figure;
L = linkage( pdist(pdf_sp, @(Xi, Xj)Hhellinger(Xi, Xj)), 'weighted' );
dendrogram(L, 79, 'Labels', folders5, 'Orientation', 'right', 'ColorThreshold', 0.45*max(L(:, 3)));
title('Phylloacoustic tree using KDE+Hellinger distance');