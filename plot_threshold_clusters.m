function [ clusters ] = plot_threshold_clusters( L, max_d, n_points, plot )
%PLOT_THRESHOLD_CLUSTERS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    plot = false;
end

mesh = 0:(1/n_points):max_d;
distances = L(:, 3);
n_init = length( distances );
clusters = zeros( n_init, 1 );
for i=1:length(mesh)
    clusters(i) = n_init - length( distances( distances < mesh(i) ) );
end

if plot
    figure;
    plot( mesh, clusters );
end

end

