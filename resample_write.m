function [  ] = resample_write( filename, ofs, fs, p, q, sections, seconds, path, name )
%RESAMPLE_WRITE Summary of this function goes here
%   Detailed explanation goes here

[ signal ofs ] = audioread( filename, [1 ofs*sections*seconds] );
signal = signal(:, 1);
x = resample( signal, p, q );
finish = 0;
for i=1:sections
    begin = finish + 1;
    finish = fs*seconds*i;
    seg = x( begin:finish );
    audiowrite( strcat( path, name, sprintf('%i', i), '.wav' ), seg, fs );
end

end

