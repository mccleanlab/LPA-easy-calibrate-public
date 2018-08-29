% SUMMARY:
% FILL IN
%
% INPUTS:
%  none
%
%
% OUTPUTS:
%  TBD
%
%
%    Written by Kieran Sweeney
%               University of Wisconsin-Madison
%               Department of Biomedical Engineering
%               1550 Engineering Drive ECB 3156
%               Madison, WI 53705
%               ksweeney2@wisc.edu
%
%
%    Last revised on August 23, 2018


%% Prepare to run script
clear all; close all; clc;

%% Set segmentation parameters
ampthresh = 0.7; % Percent of max intensity threshold for segmenting wells
sdthresh = 1.5; % Threshold for discarding unwanted data points from wells (ie, data captured when moving sensor to well)
minPeakDist = 0; % Can optionally enforce minimum distance between well peaks to improve well identification

%% Set measurement parameters
numRows = 4; % Rows A-D of 24 well plate
numColumns = 6; % Columns 1-6 of 24 well plate
numWells = numRows*numColumns;
rowNames = ['A'; 'B'; 'C'; 'D'];
channelsPerWell = 2;
cmap = lines(numWells);

%% Set calibration round
UserAnswerRound = inputdlg('Enter calibration round');
calibrationRound = str2num(UserAnswerRound{1});

%% Set calibration property
UserAnswerProperty = questdlg('Select property to calibrate','Calibration property','gcal','dc','gcal');
switch UserAnswerProperty
    case 'gcal'
        calFile = 'gcal';
        maxCal = 255; % Initial max gcal value
    case 'dc'
        calFile = 'dc';
        maxCal = 63; % Initial max dc value
end

%% Load data
[file, folder] = uigetfile('*',['Select round ' num2str(calibrationRound) ' channel 1 data']);
file_ch1 = [folder file];
[file, folder] = uigetfile('*',['Select round ' num2str(calibrationRound) ' channel 1 data']);
file_ch2 = [folder file];
files = {file_ch1, file_ch2};

[filepath,name,ext] = fileparts(files{1});

switch ext
    case '.csv'
        dataStartLine = 16;
        dataColumn = 4;
        reverseData = 0;
    case '.txt'
        dataStartLine = 2;
        dataColumn = 2;
        reverseData = 1;
end

%% Select output folder
outputFolder = uigetdir('','Select folder to store calibration files');

%% ID wells
figure('Name', ['Round ' num2str(calibrationRound) ' Well Identification']);

for c = 1:channelsPerWell
    %% Select channel data
    data = files{c};
    opts = detectImportOptions(data);
    opts.DataLine = dataStartLine;
    data = readtable(data,opts);
    data = data{:,dataColumn};
    
    if reverseData~=0
        data = wrev(data);
    end
        
    %% Preprocess data
    subplot(2,1,c); hold on;
    title(['Channel ' num2str(c)]);
    
    time = 1:length(data);
    plot(time,data(:));
    wellIntensity = data(:);
    dark = median(wellIntensity(wellIntensity<ampthresh*max(wellIntensity))); % Calculate background dark intensity
    wellIntensity = wellIntensity - dark; % Subtract out background intensity
    plot(time,wellIntensity)
    wellIntensity(wellIntensity<ampthresh*max(wellIntensity)) = nan; % Exclude intensity data from outside wells
        
    %% ID wells
    wellIntensityMask = zeros(length(wellIntensity),1);
    wellIntensityMask(wellIntensity>0) = 1;
    
    [pks, locs, width] = findpeaks(wellIntensityMask,'MinPeakDistance',minPeakDist);
    locs = locs + round(width/2);
    [val,idx] = min(abs(time-locs));
    wellID = time(idx)';
    plot(locs, wellIntensity(locs),'k*','MarkerSize',10);

    if max(idx) > numWells
        disp(['Warning! More than ' num2str(numWells) ' wells detected in channel ' num2str(c) ' row ' rowNames(i)]);
    elseif max(idx) < numWells
        disp(['Warning! Fewer than ' num2str(numWells) ' wells detected in channel ' num2str(c) ' row ' rowNames(i)]);
    end
    
    %% Calculate average well intensities and exlude outliers (captured when moving sensor to well)
    for j = 1:numWells
        m = nanmedian(wellIntensity(wellID==j));
        sd = nanstd(wellIntensity(wellID==j));
        wellIntensity(abs(wellIntensity - m)>sdthresh*sd & wellID==j) = nan;
        wellMean(c,j) = nanmean(wellIntensity(wellID==j));
        wellSD(c,j) = nanstd(wellIntensity(wellID==j));
        text(locs(j), 1.065*max(wellIntensity(:)), sprintf('%s',['Well ' num2str(j)]), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FontSize',8,'Color',cmap(j,:));
        plot(time(wellID==j), wellIntensity(wellID==j),'o','Color',cmap(j,:))
    end
    
    set(gca,'Ylim', [0, 1.2*max(wellIntensity(:))]);
end

%% Format intensity data for calibration
rawIntensity = zeros(numRows, numColumns*channelsPerWell);
rawIntensity(:,1:2:end-1) = vec2mat(wellMean(1,:),numColumns);
rawIntensity(:,2:2:end) = vec2mat(wellMean(2,:),numColumns);

%% Refine calibration values from last round (if applicable)
if calibrationRound > 1
    calPrevious = csvread([strtrim(outputFolder) '\' strtrim(calFile) '_round_' num2str(calibrationRound - 1) '.csv']);
    intensities = rawIntensity./(calPrevious/maxCal);
else
    intensities = rawIntensity;
end

minIntensity = min(intensities(:));
relIntensity = intensities/minIntensity;
cal = 1./relIntensity;
cal = round(cal*maxCal);

%% Save calibration values
dlmwrite([outputFolder '\' calFile '_round_' num2str(calibrationRound) '.csv'],cal, 'delimiter', ',', 'precision', 9);
dlmwrite([outputFolder '\rawIntensities_round_' num2str(calibrationRound) '.csv'],rawIntensity, 'delimiter', ',', 'precision', 9);

%% Calculate results
relMaxIntensity = intensities/max(intensities(:));
relCalIntensity = intensities.*cal/max(max(intensities.*cal));
plateMean = mean(rawIntensity(:));
plateSD = std(rawIntensity(:));
plateCV = plateSD/plateMean;

%% Display and plot results
disp(['Round ' num2str(calibrationRound) ' mean plate intensity = ' num2str(plateMean*1E6) ' uW']);
disp(['Round ' num2str(calibrationRound) ' SD plate intensity = ' num2str(plateSD*1E6) ' uW']);
disp(['Round ' num2str(calibrationRound) ' CV = ' num2str(plateCV*100) '%']);

figure('Name', ['Round ' num2str(calibrationRound) ' calibration'])
plateData = {rawIntensity*1E6, cal};
titles = {['Round ' num2str(calibrationRound) ' raw intensities (uW)'], ['Round ' num2str(calibrationRound) ' ' calFile ' values']};

for i = 1:2
    subplot(2,1,i)
    heatmap(plateData{i});
    colorbar;
    title(titles{i});
end

clearvars