init;

n = randi(length(folders));
s = strcat('../../dataset/', folders{n}, '/*.wav');
figure;
hold on;
i=1;

for file=dir( s )'
    f{i} = file.name;
    kde( formants( strcat('dataset/', folders{n}, '/', file.name) ) );
    i = i + 1;
end

legend(f);
clear f;