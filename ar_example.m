function [ params ] = ar_features( folder, p )
    %Creates a p-th order AR model
    model = arima( 'ARLags', 1:p, 'Constant', 0 );

    %Reads a signal and estimates the parameters for its p-th order AR model
    files = dir( strcat( folder, '/*.wav' ) );
    i = 1;
    params = zeros( length(files), p );

    for file = files'
        [ y Fs ] = audioread( strcat( folder, '/', file.name ) );
        EstMdl = estimate( model, y, 'Display', 'off' );
        params(i, :) = cell2mat( EstMdl.AR );
        i = i+1;
    end
end