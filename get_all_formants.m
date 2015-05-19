load folders

for i=1:length(folders)
    F{i} = folder_formants(strcat('../../dataset/', folders{i}))
end

mapObj = containers.Map(folders',F);