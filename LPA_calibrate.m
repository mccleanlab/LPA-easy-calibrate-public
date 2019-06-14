% NOTE:
% This script is based on calibration.m by Sebastian Castillo-Hair, as
% described in Gerhardt, K. P. et al. An open-hardware platform for
% optogenetics and photobiology. Nat. Publ. Gr. 1–13 (2016).
% doi:10.1038/srep35363
%
% SUMMARY: This script imports files containing light intensity
% measurements acquired by a ThorLabs optical power meter, maps the
% intensities to the wells of a light plate apparatus (LPA), and calculates
% calibration values (either dc or gcal) such that each LED in the LPA
% outputs an equal amount of light. Multiple rounds of calibration may be
% needed to achieve this. If needed, users can change how intensity
% measurements are mapped to LPA wells by changing the segmentation
% parameters.
%
% The measurement files to be imported for channels 1 and 2 (corresponding
% to the top and bottom sets of LEDs, respectively) and the current
% calibration round are specified by the user via UI prompts.
%
% Whether the script calculates dc or gcal calibration values is also
% specified by UI prompt. The mapped intensity values are saved to
% meanIntensities_round_*.csv and the calibration files are saved to either
% dc_round_*.csv or gcal_round_*.csv. The .csv files are saved to an output
% folder specified by UI prompt. A user may also choose not to export
% calibration values and simply measure the light output of the LPA. In
% this case the user can select whether to measure the top and bottom sets
% of LEDs independently (as during calibration) or together (to measure net
% well intensity) for a specified number of wells. The script also displays
% the mean, standard deviation, and coefficient of variation of light
% intensities across the plate. Please note that the calculated intensity
% and gcal values for each channel (top and bottom LEDs) are displayed in
% Figure 2 as vertically interleaved to match the layout of the LPA, but
% are saved as horizontally interleaved values to match the formating used
% by the LPA.
%
% This script has been tested with 1) .csv files containing light intensity
% measurements acquired via ThorLabs Optical Power Monitor v1.0.2149.55
% software and 2) .txt files containing light intensity measurements
% acquired via ThorLabs PM100 Utility Version 3.0. It may not work directly
% with other file types.

%
% % INPUTS:
% calibrationRound: controls file import/export, entered via UI prompt
% Measurement files for channels 1 and 2: entered via UI prompt
%
% OUTPUTS:
% meanIntensities_round_*.csv: Intensity measurements mapped to LPA
% gcal_round_*.csv: LPA gcal values
% dc_round_*.csv: LPA dc values
%
% WRITTEN BY:
%   Kieran Sweeney
%   University of Wisconsin-Madison
%   Department of Biomedical Engineering
%   1550 Engineering Drive ECB 3156
%   Madison, WI 53705
%   ksweeney2@wisc.edu
%
% Last revised on February 21, 2019

%% Prepare to run script
clearvars ; close all; clc;

%% Set segmentation parameters
ampthresh = 0.5; % Fraction of max intensity threshold for segmenting wells
% sdthresh = 2; % Threshold for discarding unwanted data points from wells (ie, data captured when moving sensor to well)
minPeakDist = 0; % Can optionally enforce minimum distance between well peaks to improve well identification

%% Measurement parameters
numRows = 4; % Rows of 24 well plate
rowNames = {'A' 'B' 'C' 'D'};
rowLabels = {'A_top', 'A_bottom','B_top', 'B_bottom','C_top', 'C_bottom','D_top', 'D_bottom'};
numColumns = 6; % Columns of 24 well plate
columnNames = {'1' '2' '3' '4' '5' '6'};
numWells = numRows*numColumns;
sensorUnits = 'uW/cm^2';

%% Set calibration property by UI prompt
calibrationType = questdlg('Select property to calibrate','Calibration type','gcal','dc','none','none');
switch calibrationType
    case 'gcal'
        calFile = 'gcal';
        maxCal = 255; % Initial max gcal value
        UserAnswerRound = inputdlg('Enter calibration round');
        calibrationRound = str2double(UserAnswerRound{1});
        importPrompt = ['Select round ' num2str(calibrationRound)];
        indepLEDs = 'Yes';
    case 'dc'
        calFile = 'dc';
        maxCal = 63; % Initial max dc value
        UserAnswerRound = inputdlg('Enter calibration round');
        calibrationRound = str2double(UserAnswerRound{1});
        importPrompt = ['Select round ' num2str(calibrationRound)];
        indepLEDs = 'Yes';
    case 'none'
        indepLEDs = questdlg('Measure top and bottom LEDs independently?','Measurement type','Yes','No','No');
        numWells = inputdlg('Enter number of wells measured');
        numWells = str2double(numWells{1});
        if numWells~=24
            numRows = numWells;
            numColumns = 1;
        end
        measureOnly = 1;
        calibrationRound = 0;
        importPrompt = 'Select';
end

%% Import intensity measurements for channels 1 and 2 by UI prompt
switch indepLEDs
    case 'Yes'
        channelsPerWell = 2;
        [file, folder] = uigetfile('*',[importPrompt ' channel 1 data']);
        file_ch1 = [folder file];
        [file, folder] = uigetfile('*',[importPrompt ' channel 2 data']);
        file_ch2 = [folder file];
        files = {file_ch1, file_ch2};
    case 'No'
        channelsPerWell = 1;
        [file, folder] = uigetfile('*',[importPrompt ' data']);
        files = {[folder file]};
end

% Interpret measurements based on file extension
[filepath,name,ext] = fileparts(files{1});

switch ext
    case '.csv' % Finds intensity measurements for files from ThorLabs Optical Power Monitor v1.0.2149.55
        dataStartLine = 16; % Line where intensity measurements start
        dataColumn = 4; % Column of imported data containing intensity measurements
    case '.txt' % Finds intensity measurements in files from ThorLabs PM100 Utility Version 3.0
        dataStartLine = 2; % Line where intensity measurements start
        dataColumn = 2; % Column of imported data containing intensity measurements
        reverseData = 1; % Reverses order of measurements from ThorLabs PM100 Utility Version 3.0
end

%% Select output folder by UI prompt
outputFolder = uigetdir('','Select output folder for measurement and calibration files');

%% Identify wells and measure light intensity per well
figure('Name', ['Round ' num2str(calibrationRound) ' Well Identification']);
cmap = lines(numWells);
wellMean = nan(channelsPerWell,numWells);
wellSD = nan(channelsPerWell,numWells);

for c = 1:channelsPerWell
    % Create subplot per channel
    subplot(2,1,c); hold on;
    title(['Channel ' num2str(c)],'Interpreter', 'none');
    
    % Extract measurement data from file
    warning('off','MATLAB:textio:io:UnableToGuessFormat');
    data = files{c};
    opts = detectImportOptions(data);
    opts.DataLine = dataStartLine;
    data = readtable(data,opts);
    data = data{:,dataColumn};
    
    % Reverse measurement data if necessary
    if exist('reverseData')~=0
        data = wrev(data);
    end
    
    % Preprocess measurement data by thresholding and background subtraction
    time = 1:length(data);
    plot(time,data(:));
    wellIntensity = data(:);
    dark = median(wellIntensity(wellIntensity<ampthresh*max(wellIntensity))); % Calculate background intensity
    wellIntensity = wellIntensity - dark; % Subtract out background intensity
    plot(time,wellIntensity)
    wellIntensity(wellIntensity<ampthresh*max(wellIntensity)) = nan; % Threshold to exclude intensity data from outside wells
    
    % Create binary mask signal from preprocessed measurements
    wellIntensityMask = zeros(length(wellIntensity),1);
    wellIntensityMask(wellIntensity>0) = 1;
    
    % Find and count peaks in masked signal to identify wells
    [pks, locs, width] = findpeaks(wellIntensityMask,'MinPeakDistance',minPeakDist); % Find peaks from mask
    locs = locs + round(width/2); % Center peaks
    [val,idx] = min(abs(time-locs));
    wellID = time(idx)';
    plot(locs, wellIntensity(locs),'k*','MarkerSize',10);
    
    if max(idx) > numWells
        disp(['Warning! More than ' num2str(numWells) ' wells detected in channel ' num2str(c) newline 'Try changing the ampthresh paramater']);
    elseif max(idx) < numWells
        disp(['Warning! Fewer than ' num2str(numWells) ' wells detected in channel ' num2str(c) newline 'Try changing the ampthresh paramater']);
    end
    
    % Calculate intensities for each identified well
    for j = 1:numWells
        % Exclude outliers (acquired when moving sensor to well)
        wellIntensity0 = [];
        wellIntensity0 = wellIntensity(wellID==j);
        outliers = wellIntensity0(isoutlier(wellIntensity0));
        wellIntensity(ismember(wellIntensity,outliers) & wellID==j) = nan;
        
        % Calculate well intensities
        
        wellMean(c,j) = nanmean(wellIntensity(wellID==j));
        wellSD(c,j) = nanstd(wellIntensity(wellID==j));
        count = wellIntensity(wellID==j);
        count = length(count(~isnan(count)));
        wellCount(c,j) = count;
        if numRows==4 && numColumns==6
            text(locs(j), 1.065*max(wellIntensity(:)), sprintf('%s',[rowNames{floor((j-1)/numColumns) + 1} num2str(rem(j+numColumns-1,numColumns) + 1)]), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FontSize',8,'Color',cmap(j,:));
        else
            text(locs(j), 1.065*max(wellIntensity(:)), sprintf('%s',[num2str(j)]), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FontSize',8,'Color',cmap(j,:));
        end
        
        plot(time(wellID==j), wellIntensity(wellID==j),'o','Color',cmap(j,:)); hold on
        ylabel ('Irradiance W/cm^2')
        xlabel ('Samples')
    end
    
    set(gca,'Ylim', [0, 1.2*max(wellIntensity(:))]);
end

%% Map measurements to LPA wells
% Convert well intensities vector into matrix corresponding to LPA location
meanIntensities = zeros(numRows, numColumns*channelsPerWell);
meanIntensitiesDisplay = zeros(numRows*channelsPerWell, numColumns);

if channelsPerWell == 1
    meanIntensities = vec2mat(wellMean(:),numColumns);
    meanIntensitiesDisplay = meanIntensities;
else
    meanIntensities(:,1:2:end-1) = vec2mat(wellMean(1,:),numColumns);
    meanIntensities(:,2:2:end) = vec2mat(wellMean(2,:),numColumns);
    meanIntensitiesDisplay(1:2:end-1,:) = vec2mat(wellMean(1,:),numColumns);
    meanIntensitiesDisplay(2:2:end,:) = vec2mat(wellMean(2,:),numColumns);
end

%% Calculate and display plate statistics
plateMean = mean(meanIntensities(:));
plateSD = std(meanIntensities(:));
plateCV = plateSD/plateMean;

disp(['Round = ' num2str(calibrationRound)])
disp(['Mean LPA intensity = ' num2str(plateMean*1E6) ' ' sensorUnits]);
disp(['SD LPA intensity = ' num2str(plateSD*1E6) ' ' sensorUnits]);
disp(['LPA CV = ' num2str(plateCV*100) '%']);

%% Export measurements only (if 'none' selected in calibration type prompt)
if exist('measureOnly')~=0
    dlmwrite([outputFolder '\meanIntensities_round_' num2str(calibrationRound) '.csv'],meanIntensities, 'delimiter', ',', 'precision', 9);
    % Plate mapped measurements
    figure('Name', ['Measured intensities'])
    
    if channelsPerWell==1
        rowLabels = rowNames;
    end
    
    if numWells~=24
        h = heatmap(meanIntensitiesDisplay*1E6);
    else
        h = heatmap(columnNames, rowLabels, meanIntensitiesDisplay*1E6);
    end
    
    colorbar;
    title(['Mean intensities (' sensorUnits ')'])
    %     clearvars -except rmeanIntensities cal plateMean plateSD plateCV
    return  % Stop script here
end

%% Calculate calibration values

% Scale data with information from previous calibration rounds (if applicable)
if calibrationRound > 1
    % Set heatmap colorbar range to span initial measured intensities
    intensitiesInit = csvread([strtrim(outputFolder) '\meanIntensities_round_1.csv']);
    colorbarMin = min(intensitiesInit(:));
    colorbarMax = max(intensitiesInit(:));
    % Scale intensities by calibration values from previous round (if applicable)
    calPrevious = csvread([strtrim(outputFolder) '\' strtrim(calFile) '_round_' num2str(calibrationRound - 1) '.csv']);
    intensities = meanIntensities./(calPrevious/maxCal);
else
    intensities = meanIntensities;
end

% Calculate calibration values
minIntensity = min(intensities(:));
relIntensity = intensities/minIntensity;
cal = 1./relIntensity;
cal = round(cal*maxCal);

% Reformat calibration values for display
calDisplay = zeros(8,6);
calDisplay(1:2:end-1,:) = cal(:,1:2:end-1);
calDisplay(2:2:end,:) = cal(:,2:2:end);

%% Output results
% Save mapped measurements and calibration values
dlmwrite([outputFolder '\' calFile '_round_' num2str(calibrationRound) '.csv'],cal, 'delimiter', ',', 'precision', 9);
dlmwrite([outputFolder '\meanIntensities_round_' num2str(calibrationRound) '.csv'],meanIntensities, 'delimiter', ',', 'precision', 9);

%% Plot results
figure('Name', ['Round ' num2str(calibrationRound) ' calibration'])

% Plate mapped measurements
subplot(2,1,1);
h1 = heatmap(columnNames, rowLabels, meanIntensitiesDisplay*1E6);
colorbar;
if calibrationRound > 1
    h1.ColorLimits = [colorbarMin colorbarMax]*1E6;
end
title(['Round ' num2str(calibrationRound) ' mean intensities (' sensorUnits ')'])

% Plot calibration values
subplot(2,1,2);
h2 = heatmap(columnNames, rowLabels, calDisplay);
colorbar;
title(['Round ' num2str(calibrationRound) ' ' calFile ' values'])

clearvars -except meanIntensities cal plateMean plateSD plateCV