clear all; close all; clc;

%% Set calibration round and channel
calibrationRound = 4;

%% Select source data folders
folder_ch1 = 'D:\Kieran\MATLAB_Scripts\LPA\LPA04\rawData_calibrate\Round4_channel1\';
folder_ch2 = 'D:\Kieran\MATLAB_Scripts\LPA\LPA04\rawData_calibrate\Round4_channel2\';

%% Select output folder
outputFolder = 'D:\Kieran\MATLAB_Scripts\LPA\LPA04\gcal';

%% Set parameters
ampthresh = 0.7; % Percent of max intensity threshold for segmenting wells
sdthresh = 1.5; % Threshold for discarding unwanted data points from wells (ie, data captured when moving sensor to well)
minPeakDist = 0;
maxCal = 255; % Initial max value of LED intensity

numRows = 4; % Rows A-D of 24 well plate
numColumns = 6; % Columns 1-6 of 24 well plate
rowNames = ['A'; 'B'; 'C'; 'D'];
channelsPerWell = 2;
totalColumns = numColumns*channelsPerWell;

%% Measure intensity from raw power meter data and assign to well
folders = {folder_ch1, folder_ch2};

for c = 1:channelsPerWell
    folder = folders{c};
    figure('Name', ['Channel ' num2str(c) ' well identification']);
    %% Import data from Thorlabs power meter
    A = readtable([folder 'A.txt']); A = wrev(A{:,2});
    B = readtable([folder 'B.txt']); B = wrev(B{:,2});
    C = readtable([folder 'C.txt']); C = wrev(C{:,2});
    D = readtable([folder 'D.txt']); D = wrev(D{:,2});
    
    data = nan(100,4);
    data(1:length(A),1) = A;
    data(1:length(B),2) = B;
    data(1:length(C),3) = C;
    data(1:length(D),4) = D;
    
    %% Measure light intensity per well
    for i = 1:numRows
        time = 1:length(data(:,i));
        subplot(4,1,i); hold on;
        plot(time,data(:,i)); title(['Row ' rowNames(i) ' channel ' num2str(c)]); ylabel('Intensity (W)')
        wellIntensity = data(:,i);
        dark = median(wellIntensity(wellIntensity<ampthresh*max(wellIntensity))); % Calculate background dark intensity
        wellIntensity = wellIntensity - dark; % Subtract out background intensity
        plot(time,wellIntensity)
        wellIntensity(wellIntensity<ampthresh*max(wellIntensity)) = nan; % Exclude intensity data from outside wells
        
        
        %% ID wells
        wellIntensityMask = zeros(length(wellIntensity),1);
        wellIntensityMask(wellIntensity>0)=1;
        
        [pks, locs, width] = findpeaks(wellIntensityMask,'MinPeakDistance',minPeakDist);
        locs = locs + round(width/2);
        [val,idx] = min(abs(time-locs));
        wellID = time(idx)';
        plot(locs, wellIntensity(locs),'g*');
        plot(time,(wellID/max(wellID))*max(wellIntensity));
        if max(idx) ~= numColumns
            disp(['Warning! More than ' num2str(numColumns) ' wells detected in channel ' num2str(c) ' row ' rowNames(i)]);
        end
        
        %% Calculate average well intensities and exlude low outliers (captured when moving sensor to well)
        for j = 1:numColumns
            m = nanmedian(wellIntensity(wellID==j));
            sd = nanstd(wellIntensity(wellID==j));
            wellIntensity(wellIntensity<(m-sdthresh*sd) & wellID==j) = nan;
            wellMean(c,i,j) = nanmean(wellIntensity(wellID==j));
            wellSD(c,i,j) = nanstd(wellIntensity(wellID==j));
        end
        
        plot(time,wellIntensity,'ro');
        
    end
end

%% Format measured intensity data for calibration
rawIntensity = zeros(numRows, totalColumns);
rawIntensity(:,1:2:end-1) = wellMean(1,:,:);
rawIntensity(:,2:2:end) = wellMean(2,:,:);

%% Refine values from last round of calibration (if applicable)
if calibrationRound > 1
    calPrevious = csvread([strtrim(outputFolder) '\gcal_round_' num2str(calibrationRound - 1) '.csv']);
    intensities = rawIntensity./(calPrevious/maxCal);
else
    intensities = rawIntensity;
end

minIntensity = min(intensities(:));
relIntensity = intensities/minIntensity;
cal = 1./relIntensity;
cal = round(cal*maxCal);

%% Save calibration values
dlmwrite([outputFolder '\gcal_round_' num2str(calibrationRound) '.csv'],cal, 'delimiter', ',', 'precision', 9);

%% Plot results

relMaxIntensity = intensities/max(intensities(:));
relCalIntensity = intensities.*cal/max(max(intensities.*cal));
plateMean = mean(rawIntensity(:));
plateSD = std(rawIntensity(:));
plateCV = plateSD/plateMean;

disp(['Round ' num2str(calibrationRound) ' mean plate intensity = ' num2str(plateMean*1E6) ' uW']);
disp(['Round ' num2str(calibrationRound) ' SD plate intensity = ' num2str(plateSD*1E6) ' uW']);
disp(['Round ' num2str(calibrationRound) ' CV = ' num2str(plateCV*100) '%']);

figure('Name', ['Round ' num2str(calibrationRound) ' calibration'])
data = {rawIntensity*1E6, cal};
titles = {['Round ' num2str(calibrationRound) ' raw intensities (uW)'], ['Round ' num2str(calibrationRound) ' calibration values']};

for i = 1:2
    subplot(2,1,i)
    heatmap(data{i});
    colorbar;
    title(titles{i});
end

clearvars -except calibrationRound outputFolder maxCal numRows totalColumns wellMean wellSD