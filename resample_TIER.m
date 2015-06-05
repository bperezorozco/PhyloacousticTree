load folders;
tfs = 44100;

for folder=folders'
    s = strcat('../../dataset/', folder{1}, '/*.wav');
    
    for file=dir( s )'
        filename = strcat('../../dataset/', folder{1}, '/', file.name);
        [x fs] = audioread( filename );
        
        if fs ~= tfs
            fs
        end
        %y = resample(x, fs, tfs);
        %audiowrite( filename, y, tfs );
    end

end