%% importPIVdata_ASCII -  Makes one .mat file out of multiple ASCII files.
%%                        The script must be located in the Codes folder to
%%                        properly work.
%% ------------------------------------------------------------------------
%  inputFolder- Name of the directory where the ASCII files containing data for
%           each timestep is located. The mat file is stored where this
%           folder is located. Also the ascii files must be in a folder
%           called "Export_ASCII"
% zvel          -   Boolean used to determine stereo PIV
%% ========================================================================

function importPIVdata_ASCII(inputFolder, zvel)

dirname = sprintf("%s/Export_ASCII",inputFolder);

Data = struct;  %Struct to keep all data

% Vector with the names of the.txt files
txtfiles = dir(fullfile(dirname,'*.txt'));
numfiles = length(txtfiles);

%Write velcome message
fprintf("\n=============================================================");
fprintf("\nImporting PIV data ""%s""...\n\n", dirname);
% frequency = input("\nCapture frequency [Hz]: ");
frequency = 15;

%% Run through all the files and store the data in a mat file
for i = 1:1:numfiles
    % Write out the progress to the user
    progress = (i/numfiles)*100;
    fprintf("\nImporting file %i of %i (%3.0f%%)", i, numfiles, progress);
    
    filename = sprintf('%s\\%s',txtfiles(i).folder,txtfiles(i).name);
    dat = importdata(filename);
    
    % For first file determine lenght of vectors as we are assuming
    % uniformly gridded data
    if(i == 1)
        xl = checkVector(dat.data(:,1)); %Checks length of x-coordinates
        yl = length(dat.data(:,2))/xl;   %Based on xl, find length of
        %y vector)
        xvec = dat.data(:,1);            %Allocate values
        yvec = dat.data(:,2);
        %Save coordinate info to the output struct
        Data(1).freq = frequency;
        Data(1).xl = xl;
        Data(1).yl = yl;
        Data(1).x = xvec(1:xl);
        Data(1).y = yvec(1:xl:end);
    end
    
    % Get all the velocity data
    vx = dat.data(:,3);
    vy = dat.data(:,4);
    
    if zvel
        vz = dat.data(:,5);
    end
    
    % Go through veklocity data and place in matrices used for contour
    % plotting
    c = 0;
    velXArray = zeros(yl,xl);
    velYArray = zeros(yl,xl);
    if zvel
        velZArray = zeros(yl,xl);
    end
    for y = 1:yl
        for x = 1:xl
            c = c+1;
            velXArray(y,x) = vx(c);
            velYArray(y,x) = vy(c);
            if zvel
                velZArray(y,x) = vz(c);
            end
            
        end
    end
    
    % Save to  a struct
    Data(i).velXArray = velXArray;
    Data(i).velYArray = velYArray;
    if (zvel)
        Data(i).velZArray = velZArray;
    end
end

%Save to  a mat file
temp = char(inputFolder);   %To extract only 4 first characters
matFilename = sprintf("%s/%s.mat",inputFolder,temp(end-3:end));
fprintf("\n\n...Saving data to %s...", matFilename);
save(matFilename,'Data','-v7.3');

fprintf("\n\nFinished importing data!");
fprintf("\n=============================================================");
end




