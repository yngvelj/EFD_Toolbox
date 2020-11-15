%Force measurement post-processing
%Frequency analysis and Drag/Lift coefficient calculations
%%
clear all;
%%
path = uigetdir();%Choose the folder with the Kraftmålinger folder in the prosjektoppgave folder
cd(path);


matO = dir('1000/*.mat'); %Osillatory
matP = dir('7600/*.mat'); %Plate without bump
mat7000 = dir('7000/*.mat');
mat7430 = dir('7430/*.mat');
mat7440 = dir('7440/*.mat');
mat7450 = dir('7450/*.mat');
mat7460 = dir('7460/*.mat');
mat7470 = dir('7470/*.mat');
mat7480 = dir('7480/*.mat');
mat7490 = dir('7490/*.mat');
mat7500 = dir('7500/*.mat');
mat7510 = dir('7510/*.mat');

matvec.Data = {matO,matP,mat7000,mat7430,mat7440,mat7450,mat7460,mat7470,mat7480,mat7490,mat7500,mat7510};
counter = 0;
for i = 1:12
    for j = 1:length(matvec.Data{i})
        counter = counter +1;
        run = matvec.Data{i}(1).name(1:4);
        runnr = matvec.Data{i}(j).name(1:4);
        Data(counter).run = run;
        Data(counter).name = runnr;
        Data(counter).data = load(sprintf('%s/%s',matvec.Data{i}(j).folder,matvec.Data{i}(j).name));
    end
    
end
save("ForceData.mat",'Data','-v7.3');

%%
%LDV-data:
ldvpath = uigetdir(); %Choose  ldv-Data folder in prosjektoppgave
cd(ldvpath);
folders = dir();
counter = 0;
for i = 1:length(folders)
    if(all(ismember(folders(i).name,'0123456789')))%Check if folder acutally is a run-number
        counter = counter +1;
        LDVdata(counter).name = folders(i).name;
        samplesname = sprintf('%s\\%s\\%s.SPEED.MSEBP.txt',folders(i).folder,folders(i).name,folders(i).name);
        LDVdata(counter).samples = importLDVsamples(samplesname);
        statname = sprintf('%s\\%s\\%s.STATISTICS.MSEBP.txt',folders(i).folder,folders(i).name,folders(i).name);
        LDVdata(counter).stats = ImportLDVStat(statname);
    end
    
end
%%
save("LDVdata.mat",'LDVdata','-v7.3');




%% Save data to a mat-file
% save("Data.mat",'Data','-v7.3');



%Function for extracting LDV-time samples from miniLDV(Freestream
%measurements)
function LDVSamples = importLDVsamples(filename)

Data = fileread(filename); %Read data from file

Data = strrep(Data, ',', '.'); %Change from komma to dot
fid = fopen('data.txt', 'w+');
fwrite(fid, Data, 'char'); %skriver dataen til en fil
A = importdata('data.txt'); %loads the data into one variabel

LDVSamples.timevec = A(:,1);
LDVSamples.Speedvec = A(:,3);
LDVSamples.SNRvec = A(:,4);


end


%Function for importing statistics file LDV(freestream) miniLDV
function LDVStats = ImportLDVStat(filename)
Data = fileread(filename); %Read data from file

Data = strrep(Data, ',', '.'); %Change from komma to dot
fid = fopen('data.txt', 'w+');
fwrite(fid, Data, 'char'); %skriver dataen til en fil
A = importdata('data.txt'); %loads the data into one variabel

LDVStats.mean =A(1);
LDVStats.Speed_rms = A(2);
LDVStats.SNR = A(3);
LDVStats.Nsamples = A(6);
LDVStats.duration = A(7);

end


