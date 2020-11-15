%% Fishing net PIV
%% ========================================================================
% Project: Fishing net in CWT
% Date: Summer 2020
% Subject: streamwise mono PIV to investigate jets btwn twines
% Relevant runs: 1000 series, 3000 series
%% ========================================================================
% Remember to first create a .mat file using "importPIVdata_ASCII"!
clear all;
close all;
clc;
%% Import data
% Create  a data struct with relevant PIV and LDV data
% Input:
% -------------------------------------------------------------------------
runSeries = "3000";     %The general series for the runs (1000)
runNumber = "3110";     %Run number
viscosity = 1.1175e-6;       %Viscosity of the run
d = 3;                  %Diameter of netting twine [mm]
zvel = 0;               %1 for stereo PIV (0 for 1000-series)
timestep = 1;           %The timestep to look at
% -------------------------------------------------------------------------

% Find the folders where things are
currentFolder = pwd;                        %Full file path of current folder
idcs = strfind(currentFolder,filesep);      %indexes of file seperators
aboveFolder = currentFolder(1:idcs(end)-1); %Above folder

%Figure folder for run
figureFolder = sprintf("%s\\%s\\Figures",aboveFolder,runSeries);
%Check if folder exists and create new folder if it doesnt
if ~isfolder(figureFolder)
    newFolder = sprintf("%s\\%s",aboveFolder,runSeries);
    mkdir(newFolder,"Figures");
end

%piv file: import from .mat file into workspace for faster plotting
filename_piv = sprintf("%s\\%s\\PIV\\%s\\%s.mat",aboveFolder,runSeries,...
    runNumber,runNumber);
data_piv = importdata(filename_piv);

%ldv file: Import avg properties from the run from the STATISTICS file
filename_ldv = sprintf("%s\\%s\\LDV\\%s\\%s",aboveFolder,runSeries,runNumber,...
    runNumber);
data_ldv = importLDVdata_SinglePointAVG(filename_ldv);
U = data_ldv.meanVel;

%% Calculate other relevant properties from the data
net = 68.1;           %Position of net relative to the PIV plane. This is
%found in the log for zero velocity
zeroArea = 0;    %Distance not part of PIV area that Davis for some
%reason have added to the edges

% Shift the x-position so its relative to the position of the net
shift = -data_piv(1).x(1)-net-zeroArea; %Amount to shift x-coordinates
data_piv(1).x = data_piv(1).x+shift;    %Shift the x-coordinates

% Calculate reynolds number
Rn = (data_ldv.meanVel*d*1e-3)./viscosity;

% Calculate quantities for each timestep
for i = 1:length(data_piv)
    % Absolute velocity
    data_piv(i).velAbsArray = ...
        sqrt(data_piv(i).velXArray.^2+data_piv(i).velYArray.^2);
    if zvel
        data_piv(i).velAbsArray = sqrt(data_piv(i).velXArray.^2 + ...
            data_piv(i).velYArray.^2+data_piv(i).velZArray.^2);
        
        % Non-dimensionalize z-velocities
        data_piv(i).velZArray = data_piv(i).velZArray./data_ldv.meanVel;
    end
    
    %Time
    data_piv(i).time = i/data_piv(1).freq;
    
end

%% Contour of absolut velocity
figure(1);
contourf(data_piv(1).x./d,data_piv(1).y./d,...
    data_piv(timestep).velAbsArray./U);
strRn = sprintf("Rn = %.0f",Rn);
title(["Absolute velocity behind net",strRn]);
xlabel("Distance behind net {x}/{d_{twine}}");
ylabel("y-position {y}/{d_{twine}}");
% Add a textbox with current timestep
strTime = sprintf("t = %.4f s",data_piv(timestep).time);
annotation("textbox", [0.6,0.1 0.1,0.1],"String",strTime,...
    "FitBoxToText","on", "BackgroundColor", "white");
% Add a colorbar
c = colorbar;
c.Label.String = '{U}/{U_\infty}';

% Save pictures to file
filename1 = sprintf("%s\\%s_situation_velocityAbs",figureFolder,runNumber);
saveas(1,filename1,"epsc");
saveas(1,filename1,"fig");
%% Jet Probes
% Probe positions
probePos = [-42.8;-23.6;-5.7;24.3];         %Start position of probes
ydiff = 0;                          %Difference to end position
probePos(:,2) = probePos(:,1)+ydiff;%Calculate end position based on ydiff

% plot a contour plot of the absolute velocity to plot probe positions on

figure(2);
contourf(data_piv(1).x./d,data_piv(1).y./d,...
    data_piv(timestep).velAbsArray./U);
strRn = sprintf("Rn = %.0f",Rn);
title(["Probe positions",strRn]);
xlabel("Distance behind net {x}/{d_{twine}}");
ylabel("y-position {y}/{d_{twine}}");
% Add a textbox with current timestep
strTime = sprintf("t = %.4f s",data_piv(timestep).time);
annotation("textbox", [0.6,0.1 0.1,0.1],"String",strTime,...
    "FitBoxToText","on", "BackgroundColor", "white");
% Add a colorbar
c = colorbar;
c.Label.String = "{U}/{U_\infty}";

nProbes = size(probePos,1);
probeData = cell(1,nProbes);
for i = 1 :nProbes
    fprintf("\nSampling probe %i...",i);
    fprintf("\n-------------------------------------------------------\n");
    probeData{i} = getLinearProbeData(data_piv,d,probePos(i,:),2,i);
    fprintf("\nSampling probe %i is done!",i);
    fprintf("\n-------------------------------------------------------\n");
end

legend("show");

% Save pictures to file
filename2 = sprintf("%s\\%s_situation_probePositions",figureFolder,runNumber);
saveas(2,filename2,"epsc");
saveas(2,filename2,"fig");

%% Plot the probe data

% Instantenous velocities
figure(3)
hold on;
for i = 1 :nProbes
    lstr = sprintf("Probe %i ",i);
    plot(probeData{i}(1).x./d,probeData{i}(timestep).velX./U,"DisplayName",...
        strcat(lstr, "x-velocity"));
    plot(probeData{i}(1).x./d,probeData{i}(timestep).velY./U,"DisplayName",...
        strcat(lstr, "y-velocity"));
end
hold off;

tstr = sprintf("Jet velocities \nRn = %.0f",Rn);
title(tstr);
xlabel("Distance behind net {x}/{d_{twine}}");
ylabel("{U}/{U_\infty}");
legend("show");
% Add a textbox with current timestep
strTime = sprintf("t = %.4f s",data_piv(timestep).time);
annotation("textbox", [0.6,0.1 0.1,0.1],"String",strTime,...
    "FitBoxToText","on", "BackgroundColor", "white");


% Save pictures to file
filename3 = sprintf("%s\\%s_probes_instantenousVel",figureFolder,runNumber);
saveas(3,filename3,"epsc");
saveas(3,filename3,"fig");

% Average velocities
figure(4)
hold on;
for i = 1 :nProbes
    lstr = sprintf("Probe %i ",i);
    plot(probeData{i}(1).x./d,probeData{i}(1).avgVelX./U,"DisplayName",...
        strcat(lstr, "x-velocity"));
    plot(probeData{i}(1).x./d,probeData{i}(1).avgVelY./U,"DisplayName",...
        strcat(lstr, "y-velocity"));
end
hold off;

tstr = sprintf("Averege jet velocities \n Rn = %.0f", Rn);
title(tstr);
xlabel("Distance behind net {x}/{d_{twine}}");
ylabel("{U}/{U_\infty}");
legend("show");

% Save pictures to file
filename4 = sprintf("%s\\%s_probes_averageVel",figureFolder,runNumber);
saveas(4,filename4,"epsc");
saveas(4,filename4,"fig");

% Should/Could place probes other places as well...

%% Vorticity

% Calculate vorticity field
omega = calculateVorticity(data_piv);

%% Plot vorticity
figure(5);
contourf(data_piv(1).x./d,data_piv(1).y./d,omega{timestep}.*(d*1e-3/U),15);

tstr = sprintf("Vorticity \n Rn = %.0f", Rn);
title(tstr);
xlabel("Distance behind net {x}/{d_{twine}}");
ylabel("y-position {y}/{d_{twine}}");
% Add a textbox with current timestep
strTime = sprintf("t = %.4f s",data_piv(timestep).time);
annotation("textbox", [0.6,0.1 0.1,0.1],"String",strTime,...
    "FitBoxToText","on", "BackgroundColor", "white");
c = colorbar;
c.Label.Interpreter = "latex";
c.Label.String = "$\frac{\omega_z d}{U_{\infty}}$";
c.Label.FontSize = 16;

% Save pictures to file
filename5 = sprintf("%s\\%s_situation_vorticity",figureFolder,runNumber);
saveas(5,filename5,"epsc");
saveas(5,filename5,"fig");

%% Power spectrum
powerSpectra(data_piv,d*15,d,U,6);
filename6 = sprintf("%s\\%s_spectrum",figureFolder,runNumber);
saveas(7,strcat(filename6," - X"),"fig");
saveas(8,strcat(filename6," - Y"),"fig");

%% Make a video
figure(9);
% Make a video writer
videoFilename = sprintf('%s\\%s_animation.avi',figureFolder,runNumber);
v = VideoWriter(videoFilename);
v.FrameRate = 15;
open(v);

for t=1:length(data_piv)
    fprintf("\nWriting animation to file %i",t);
    %Contour plot relevant velocity
    contourf(data_piv(1).x./d,data_piv(1).y./d,data_piv(t).velAbsArray./U,10);
    
    strRn = sprintf("Rn = %.0f",Rn);
    title(["Absolute velocity behind net",strRn]);
    xlabel("Distance behind net {x}/{d_{twine}}");
    ylabel("y-position {y}/{d_{twine}}");
    % Add a textbox with current timestep
    strTime = sprintf("t = %.4f s",data_piv(t).time);
    annotation("textbox", [0.6,0.1 0.1,0.1],"String",strTime,...
        "FitBoxToText","on", "BackgroundColor", "white");
    % Add a colorbar
    c = colorbar;
    c.Label.String = '{U}/{U_\infty}';
    
    F = getframe(9);
    writeVideo(v,F);
    
end


close(v);



