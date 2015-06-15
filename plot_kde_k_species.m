K = 5;
colours = 'brgykcp';
init;
figure;
hold on;
i=1;

for j=1:K
    n = randi(length(folders));
    s = strcat('../../dataset/', folders{n}, '/*.wav');
    folders{n}
    for file=dir( s )'
        f{i} = file.name;
        [bandwidth,density,xmesh] = kde( formants( strcat('dataset/', folders{n}, '/', file.name) ) );
        i = i + 1;
        plot(xmesh,density, strcat(colours(j), '.'));
    end
end
%legend(f);
clear f;