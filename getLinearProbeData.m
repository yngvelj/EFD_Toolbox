%% getLinearProbeData   -   Function that takes a PIV data struct
%%                          as input and returns the data along a linear
%%                          probe
%% ========================================================================
% data      -   PIV data struct as generated from importPIVdata_ASCII
% probePos  -   The start and end-pont of the linear probe in the format
%               [ystart yend]
% fig       -   The figure where the probe will be plotted for illustration
% probeData -   A datastruct with the probe data with format
%               xl:     Length of x-vector
%               x:      x-vector
%               y:      y-vector
%               velX:   x-velocity vector
%               velY:   y-velocity vector
%               velAbs: Absolut velocity vector
% -------------------------------------------------------------------------
function probeData = getLinearProbeData(data,d, probePos, fig, probenr)
%Make probePos dimensional
probePos = probePos.*d;

probeData = struct;     %initialize output
probeData(1).xl = data(1).xl;  %Length of x vector
probeData(1).x = (data(1).x)'; %Use all x-positions of original data

%Slope of the lines
a = (probePos(2)-probePos(1))/(probeData(1).x(end)-probeData(1).x(1));
dx = probeData(1).x(2)-probeData(1).x(1);

%Calculate y values based on slope and intercept
probeData(1).y(1) = probePos(1);
for i =2:probeData(1).xl
    probeData(1).y(i) = probeData(1).y(i-1)+a*dx;
end

%Find x and y indices of the probe-line in the data
tol= (abs(data(1).y(end))-abs(data(1).y(end-1)))/2; %Tolerance of position
%search is half of
%grid size
yind = zeros(1,probeData(1).xl);    %Pre-allocate
for i = 1:probeData(1).xl
    yind(i) = find(abs(data(1).y-probeData(1).y(i)) < tol);
end

n = length(data);

%Sample the probe
for i = 1:n
    fprintf("\nSampling probe %i for timestep %i (%3.0f%%)",probenr,i,...
        (i/n)*100);
    for j=1:probeData(1).xl
        for k=1:probeData(1).xl
            probeData(i).velX(k) = data(i).velXArray(yind(j),k);
            probeData(i).velY(k) = data(i).velYArray(yind(j),k);
            probeData(i).velAbs(k) = data(i).velAbsArray(yind(j),k);
        end
    end
end

%Find average velocities over all timesteps
probeData(1).avgVelX = zeros(1,probeData(1).xl);
probeData(1).avgVelY = zeros(1,probeData(1).xl);
probeData(1).avgVelAbs = zeros(1,probeData(1).xl);

for i = 1:n
    for j=1:probeData(1).xl
        probeData(1).avgVelX(j) = probeData(1).avgVelX(j) + ...
            probeData(i).velX(j);
        probeData(1).avgVelY(j) = probeData(1).avgVelY(j) + ...
            probeData(i).velY(j);
        probeData(1).avgVelAbs(j) = probeData(1).avgVelAbs(j) + ...
            probeData(i).velAbs(j);
    end
end

probeData(1).avgVelX = probeData(1).avgVelX./n;
probeData(1).avgVelY = probeData(1).avgVelY./n;
probeData(1).avgVelAbs = probeData(1).avgVelAbs./n;

%Plot the probe over the figure taken as input
lstr = sprintf("Probe %i",probenr);

figure(fig);
hold on;
plot(probeData(1).x./d,probeData(1).y./d, "LineWidth", 2,"DisplayName",lstr);
hold off;

end