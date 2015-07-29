%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\Bernardo\Desktop\metadata.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2015/07/11 23:08:29

%% Initialize variables.
filename = 'C:\Users\Bernardo\Desktop\metadata.csv';
delimiter = ',';

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
meta.english = dataArray{:, 1};
meta.latin = dataArray{:, 2};
meta.ord = dataArray{:, 3};
meta.fam = dataArray{:, 4};


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;