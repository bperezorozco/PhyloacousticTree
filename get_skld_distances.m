k=1;
for i=1:length(filenames)
    for j=(i+1):length(filenames)
        d(k) = skld( pdf{i}.density, pdf{j}.density )
        k=k+1;
    end
end