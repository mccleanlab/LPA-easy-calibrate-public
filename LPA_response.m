% SUMMARY:
% This script generates a response equation for a calibrated light plate
% apparatus (LPA) that relates light intensity output to IRIS value. The
% script will return the IRIS values needed to achieve the light outputs
% (in uW) listed in the variable targetLightOutputs.
% 
% The script imports files containing light intensity measurements acquired
% from an LPA with a ThorLabs power meter over a range of IRIS values. The
% measurements for each IRIS value should be saved into individual files
% with the IRIS value included in the filename. For example, one could
% measure an LPA at IRIS = 500 and log the measurements in IRIS500.csv. Do
% not include other numbers in the filename The files for each IRIS value
% can be loaded simultaneously by UI prompt. It is not necessary to measure
% the whole plate for each IRIS value (we measured 8 randomly sampled wells
% specified by the script LPA_randomizer.m), though the number of wells
% measured per IRIS value should be consistent and is specified by UI
% prompt.
% 
% This script has been tested with 1) .csv files containing light intensity
% measurements acquired via ThorLabs Optical Power Monitor v1.0.2149.55
% software and 2) .txt files containing light intensity measurements
% acquired via ThorLabs PM100 Utility Version 3.0. It may not work directly
% with other file types.
%  
% INPUTS: 
% Files containing intensity measurements for different IRIS values
%  
% OUTPUTS:
% doseTable: lists IRIS values for user-specified light intensity outputs
% doseEqn: equation relating IRIS value to light intensity output
%
% WRITTEN BY:
%   Kieran Sweeney
%   University of Wisconsin-Madison
%   Department of Biomedical Engineering
%   1550 Engineering Drive ECB 3156
%   Madison, WI 53705
%   ksweeney2@wisc.edu
%
% Last revised on August 30, 2018

%% Prepare to run script
clearvars; close all; clc;

%% Set target light outputs (in uW)
targetLightOutput = [10 20 50 75 100]'; 

%% Set segmentation parameters
ampthresh = 0.7; % Fraction of max intensity threshold for segmenting wells
sdthresh = 1; % Threshold for discarding unwanted data points from wells (ie, data captured when moving sensor to well)
minPeakDist = 0; % Can optionally enforce minimum distance between well peaks to improve well identification

%% Import intensity measurements by UI prompt
[files, folder] =  uigetfile('*','MultiSelect','on');

if ischar(files)==1
    files = {files};
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

% Extract IRIS values from filenames
numFiles = length(files);
measurementNames = erase(files,ext);
measurementValues = regexprep(files, '^(\D+)(\d+)(.*)', '$2');
measurementValues = str2double(measurementValues);

%% Set number of wells measuredby UI prompt
UserAnswerRound = inputdlg('Enter number of wells measured per IRIS value');
numWellsMeasured = str2double(UserAnswerRound{1});

%% Identify wells and measure light intensity per well
figure('Name', 'Well Identification');
cmap = lines(numWellsMeasured);
wellMean = nan(numFiles,numWellsMeasured);
wellSD = nan(numFiles,numWellsMeasured);

for i = 1:numFiles
    % Create subplot per measurement file
    subplot(numFiles,1,i); hold on;
    title(measurementNames{i},'Interpreter', 'none');
        
    % Extract measurement data from file
    file = files{i};
    opts = detectImportOptions([folder file]);
    opts.DataLine = dataStartLine;
    data = readtable([folder file],opts);
    data = data{:,dataColumn};
    
    % Reverse measurement data if necessary
    if reverseData~=0
        data = wrev(data);
    end
    
    % Preprocess measurement data by thresholding and background subtraction
    time = 1:length(data);
    plot(time,data(:));
    wellIntensity = data(:);
    dark = median(wellIntensity(wellIntensity<ampthresh*max(wellIntensity))); % Calculate background dark intensity
    wellIntensity = wellIntensity - dark; % Subtract out background intensity
    plot(time,wellIntensity)
    wellIntensity(wellIntensity<ampthresh*max(wellIntensity)) = nan; % Threshold to exclude intensity data from outside wells
    
    % Create binary mask signal from preprocessed measurements
    wellIntensityMask = zeros(length(wellIntensity),1);
    wellIntensityMask(wellIntensity>0) = 1;
    
    % Find and count peaks in masked signal to identify wells
    [pks, locs, width] = findpeaks(wellIntensityMask,'MinPeakDistance',minPeakDist);
    locs = locs + round(width/2); % Center peaks
    [val,idx] = min(abs(time-locs));
    wellID = time(idx)';
    plot(locs, wellIntensity(locs),'k*','MarkerSize',10);
    
    if max(idx) > numWellsMeasured
        disp(['Warning! More than ' num2str(numWellsMeasured) ' wells detected']);
    elseif max(idx) < numWellsMeasured
        disp(['Warning! Fewer than ' num2str(numWellsMeasured) ' wells detected']);
    end
    
    %% Calculate average well intensities and exlude outliers (captured when moving sensor to well)
    for j = 1:numWellsMeasured
        % Exclude outliers (acquired when moving sensor to well)
        m = nanmedian(wellIntensity(wellID==j));
        sd = nanstd(wellIntensity(wellID==j));
        wellIntensity(abs(wellIntensity - m)>sdthresh*sd & wellID==j) = nan;
        
        % Calculate well intensities
        wellMean(i,j) = nanmean(wellIntensity(wellID==j));
        wellSD(i,j) = nanstd(wellIntensity(wellID==j));
        text(locs(j), 1.25*max(wellIntensity(:)), sprintf('%s',['Well ' num2str(j)]), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FontSize',8,'Color',cmap(j,:));
        plot(time(wellID==j), wellIntensity(wellID==j),'o','Color',cmap(j,:))
    end
    
    set(gca,'Ylim', [0, 1.5*max(wellIntensity(:))]);
    
end

%% Generate and display response equation
% Calculate plate statistics per IRIS value
intensity = [zeros(numWellsMeasured,1), wellMean'];
measurementValues = [0; measurementValues'];
u = mean(intensity)';
sd = std(intensity)';

% Fit measured intensity versus IRIS value
[f, g] = fit(measurementValues, u,'poly1');
x = 0:1:max(measurementValues);
yfit = f.p1*x + f.p2;
doseEqn = ['Intensity = ' num2str(f.p1) '*IRIS + ' num2str(f.p2)];
disp([doseEqn newline])

% Calculate and display IRIS needed for target light output
inputIRIS = round((targetLightOutput*1E-6 - f.p2)./f.p1);
doseTable = table(inputIRIS, targetLightOutput);
doseTable.Properties.VariableNames = {'IRIS' 'Light_uW'};
disp(doseTable)

% Plot response equation versus fitted data
figure('Name', 'LPA response')
errorbar(measurementValues,u*1E6,sd*1E6,'ro','MarkerSize',2);
hold on;
h = plot(x,yfit*1E6);
title(['LPA response' newline '(R squared = ' num2str(g.rsquare) ')']); xlabel('IRIS values (AU)'); ylabel('Light intensity (uW)')
legend(h,doseEqn,'Location','northwest')
legend('boxoff')

clearvars -except doseTable doseEqn