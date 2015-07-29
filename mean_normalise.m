function [ xp ] = mean_normalise( x )
%MEAN_NORMALISE Removes mean from the signal and scales it using its std
xp = ( x - mean(x) ) / std(x);
end

