load filenames;
load folders;
load f1;
formants = containers.Map(filenames',F1);
addpath(genpath('nethelp3_2'))
addpath(genpath('netlab3_2'))
%addpath(genpath('marbox_1.1'))
addpath(genpath('hmmbox_4_1'))
%addpath(genpath('hmmbox_3_2'))
addpath(genpath('kde'))