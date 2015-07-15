init;
load f1;
figure;
hold on;

for i=0:4
    n = randi(length(filenames));
    kde(F1{n});
    s{i+1} = filenames{n};
end
legend(s);