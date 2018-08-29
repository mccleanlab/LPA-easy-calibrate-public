%% Prepare to run script
clear all; close all; clc;

%% Set target light output
lightOutput = [10 20 50 75 100]';

%% Set file parameters
numWellsMeasured = 8;

%% Set segmentation parameters
ampthresh = 0.7; % Percent of max intensity threshold for segmenting wells
sdthresh = 1; % Threshold for discarding unwanted data points from wells (ie, data captured when moving sensor to well)
minPeakDist = 0;
cmap = lines(numWellsMeasured);

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
measurementValues = regexprep(files, '^(\D+)(\d+)(.*)', '$2');
measurementValues = str2double(measurementValues);

%% Measure light intensity per well
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
    
    if max(idx) > numWellsMeasured
        disp(['Warning! More than ' num2str(numWellsMeasured) ' wells detected']);
    elseif max(idx) < numWellsMeasured
        disp(['Warning! Fewer than ' num2str(numWellsMeasured) ' wells detected']);
    end
    
    %% Calculate average well intensities and exlude outliers (captured when moving sensor to well)
    for j = 1:numWellsMeasured
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

%% Calculate plate intensity per IRIS value
intensity = [zeros(numWellsMeasured,1), wellMean'];
measurementValues = [0; measurementValues'];
u = mean(intensity)';
sd = std(intensity)';

%% Fit data
[f, g] = fit(measurementValues, u,'poly1');
x = 0:1:max(measurementValues);
yfit = f.p1*x + f.p2;
disp(['Intensity = ' num2str(f.p1) '*IRIS + ' num2str(f.p2) newline])

%% Calculate IRIS needed for target light output
inputIRIS = round((lightOutput*1E-6 - f.p2)./f.p1);
doseTable = table(inputIRIS, lightOutput);
doseTable.Properties.VariableNames = {'IRIS' 'Light_uW'};
disp(doseTable)

%% Plot data
figure;
errorbar(measurementValues,u*1E6,sd*1E6,'ro','MarkerSize',2);
hold on;
plot(x,yfit*1E6);
title(['LPA response curve' newline '(R squared = ' num2str(g.rsquare) ')']); xlabel('IRIS values (AU)'); ylabel('Light intensity (uW)')

% clearvars -except doseTable