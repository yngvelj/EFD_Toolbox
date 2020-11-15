%% Script used to run importPIVdata_ASCII function on many runs at a time

clear all;
clc;
close all;

%--------------------------------------------------------------------------
run = "3000";
%--------------------------------------------------------------------------

currentFolder = pwd;

%Go out from Code
temp = char(currentFolder);
currentFolder = temp(1:end-5);

filename = sprintf("%s/%s/PIV",currentFolder,run);

folders = dir(filename);

 for i = 3:length(folders)
    dirname = sprintf("%s/%s",filename,folders(i).name);
    importPIVdata_ASCII(dirname);
 end
