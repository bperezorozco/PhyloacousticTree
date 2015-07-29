figure;
hold on;

for i=1:5
    n = randi(length(id_folders));
    kde(F{n});
    s{i} = id_folders{n};
end
legend(s);