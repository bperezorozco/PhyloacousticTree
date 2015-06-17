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
    l = 0;
    d = 0;
    for file=dir( s )'
        f{i} = file.name;
        [bandwidth,density,xmesh] = kde( formants( strcat('dataset/', folders{n}, '/', file.name) ), 0, 8000 );
        i = i + 1;
        d = d + density;
        l = l + 1;
    end
    d = d / l;
    plot(xmesh,d, strcat(colours(j)));
    clear d;
end
%legend(f);
clear f;