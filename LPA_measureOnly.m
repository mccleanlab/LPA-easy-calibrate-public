%% Prepare to run script
clearvars; close all; clc;

%% Set file and measurement parameters
numRowsMeasured = 3;
numColumnsMeasured = 5;
numWells = numRowsMeasured*numColumnsMeasured;

%% Set segmentation parameters
ampthresh = 0.1; % Fraction of max intensity threshold for segmenting wells
sdthresh = 0.75; % Threshold for discarding unwanted data points from wells (ie, data captured when moving sensor to well)
minPeakDist = 0;
cmap = lines(numWells);

%% Load source data
[files, folder] =  uigetfile('*','MultiSelect','on');

if ischar(files)==1
    files = {files};
end

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

numFiles = length(files);
measurementNames = erase(files,ext);

%% Measure light intensity per well
figure('Name', 'Well Identification');
wellMean = nan(numFiles,numWells);
wellSD = nan(numFiles,numWells);

for i = 1:numFiles
    file = files{i};
    opts = detectImportOptions([folder file]);
    opts.DataLine = dataStartLine;
    data = readtable([folder file],opts);
    data = data{:,dataColumn};
    
    if reverseData~=0
        data = wrev(data);
    end
        
    subplot(numFiles,1,i); hold on;
    title(measurementNames{i},'Interpreter', 'none');
    
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
        disp(['Warning! More than ' num2str(numWells) ' wells detected']);
    elseif max(idx) < numWells
        disp(['Warning! Fewer than ' num2str(numWells) ' wells detected']);
    end
    
    %% Calculate average well intensities and exlude outliers (captured when moving sensor to well)
    for j = 1:numWells
        m = nanmedian(wellIntensity(wellID==j));
        sd = nanstd(wellIntensity(wellID==j));
        wellIntensity(abs(wellIntensity - m)>sdthresh*sd & wellID==j) = nan;
        wellMean(i,j) = nanmean(wellIntensity(wellID==j));
        wellSD(i,j) = nanstd(wellIntensity(wellID==j));
        text(locs(j), 1.25*max(wellIntensity(:)), sprintf('%s',['Well ' num2str(j)]), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FontSize',8,'Color',cmap(j,:));
        plot(time(wellID==j), wellIntensity(wellID==j),'o','Color',cmap(j,:))
    end
    
    set(gca,'Ylim', [0, 1.5*max(wellIntensity(:))]);
    
end

%% Display results
rawIntensitiesMean = vec2mat(wellMean,numColumnsMeasured);
rawIntensitiesSD = vec2mat(wellSD,numColumnsMeasured);
plateMean = mean(rawIntensitiesMean(:));
plateSD = std(rawIntensitiesMean(:));
plateCV = plateSD/plateMean;

disp(['Mean plate intensity = ' num2str(plateMean*1E6) ' uW']);
disp(['SD plate intensity = ' num2str(plateSD*1E6) ' uW']);
disp(['CV = ' num2str(plateCV*100) '%']);

figure('Name', 'Measured intensities');
heatmap(rawIntensitiesMean*1E6);
colorbar;
title('Measured intensities (uW)');

clearvars -except rawIntensitiesMean rawIntensitiesSD plateMean plateSD plateCV