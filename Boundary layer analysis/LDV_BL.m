%%Boundary layer plotting script
%Import LDV-data from ldv-traverse and filter each point
%The script collects velocity and position in a struct called Data
%vel, x, y, z where vel is a vector of all the samples at that point.

clear all;

path = uigetdir();%Choose the folder with LDV-station_samples ie. 7032
cd(path);

SpeedFiles = dir('Stations/*SPEED.MSEBP.txt');
StationFile = dir('*STATIONS.profile.txt');


stationStruct = readStations(StationFile.name);%Collect data about where the measurement was taken(traverse)
for i = 1:length(SpeedFiles)
    speedStruct = readSpeed(SpeedFiles(i).name);%Collect all datasamples from 

    SpeedFiles(i).name
    Data(i).speed = speedStruct.vel;%Collect each point measurement in a struct
    Data(i).x = stationStruct.x(i);
    Data(i).y = stationStruct.y(i);
    Data(i).z = stationStruct.z(i);
    Data(i).rms = stationStruct.rms(i);
    Data(i).Calculatedmean = mean(speedStruct.vel);
    Data(i).StationMean = stationStruct.meanVel(i);
end


%Plot boundaryLayer
V_inf = 0.9692;


%% 7220:
% 
BL.BL1z = [Data(1:14).z];
BL.BL1speed = [Data(1:14).Calculatedmean];

y = -flip(BL.BL1z-376.5);
u = flip(BL.BL1speed);
U = V_inf;
[delta_star, theta] = calcMomAndDisplThickness(y, u, U, 3);

BL.BL2z = [Data(31:38).z];
BL.BL2speed = [Data(31:38).Calculatedmean];

y = -flip(BL.BL2z-376.5);
u = flip(BL.BL2speed);
U = V_inf;
[delta_star, theta] = calcMomAndDisplThickness(y, u, U, 4);

figure(1)
plot(-flip(BL.BL1z-376.5),flip(BL.BL1speed./V_inf))
xlabel('y-pos[mm]');
ylabel('V/V_{inf}');
title('X = 140 mm from plate tip');

figure(2)
plot(-flip(BL.BL2z-376.5),flip(BL.BL2speed./V_inf))
xlabel('y-pos[mm]');
ylabel('V/V_{inf}');
title('X = 290 mm from plate tip');

%% 7231
% y = -flip([Data.z]-376.5);
% u = flip([Data.StationMean]);
% U = V_inf;
% [delta_star, theta] = calcMomAndDisplThickness(y, u, U, 1);
% figure(2)
% plot(-flip([Data.z]-376.5),flip([Data.StationMean]./V_inf))
% xlabel('y-pos[mm]');
% ylabel('V/V_{inf}');
% title('X = 340 mm from plate tip');

%% 7232
% y = -flip([Data.z]-376.5);
% u = flip([Data.StationMean]);
% U = V_inf;
% [delta_star, theta] = calcMomAndDisplThickness(y, u, U, 1);
% 
% figure(2)
% hold on;
% plot(-flip([Data.z]-376.5),flip([Data.StationMean]./V_inf))
% scatter(-flip([Data.z]-376.5),flip([Data.StationMean]./V_inf));
% xlabel('y-pos[mm]');
% ylabel('V/V_{inf}');
% title('X = 626 mm from plate tip');




function [speedStruct] = readSpeed(filename) 
Data = fileread(filename); %Read data from file
    if(isempty(Data))
        speedStruct.msec =  nan; %device time in milliseconds
        speedStruct.usec = nan; %device time in usec???
        speedStruct.vel = nan; %velocity
        speedStruct.SNR = nan; %signl to noise ratio
    else
        Data = strrep(Data, ',', '.'); %Change from komma to dot
        fid = fopen('data.txt', 'w+');
        fwrite(fid, Data, 'char'); %skriver dataen til en fil

        A = importdata('data.txt'); %loads the data into one variabel

        %sort the data
        speedStruct.msec = A(:,1); %device time in milliseconds
        speedStruct.usec = A(:,2); %device time in usec???
        speedStruct.vel = A(:,3); %velocity
        speedStruct.SNR = A(:,4); %signl to noise ratio
    end



end

function [stationStruct] = readStations(filename) 
    %Station definition: 1)i,   2)Axis 1 (X),   3)Axis 2 (Y),   4)Axis 3 (Z),   5)Velocity mean(m/sec),
    %   6)Velocity RMS (m/sec), 7)SNR mean, 8)SNR,  9)RMS,  10)Data rate (Hz),  
    %   11)Number acquired Acquisition length (sec),    13)Limit reached


    Data = fileread(filename); %Read data from file
    Data = strrep(Data, ',', '.'); %Change from komma to dot
    fid = fopen('stations.txt', 'w+');
    fwrite(fid, Data, 'char'); %skriver dataen til en fil

    A = importdata('stations.txt'); %loads the data into one variabel

    stationStruct.x = A(:,2);
    stationStruct.y = A(:,3);
    stationStruct.z = A(:,4);
    stationStruct.rms = A(:,6);
    stationStruct.meanVel = A(:,5);

end

