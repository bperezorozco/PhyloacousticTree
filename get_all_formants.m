load folders

for i=1:length(folders)
    H{i} = folder_formants(strcat('../../dataset/', folders{i}))
end

empObj = containers.Map(folders',H);

D_emp_2 = get_rows_distance( unfold_cell_array(values(empObj)) );