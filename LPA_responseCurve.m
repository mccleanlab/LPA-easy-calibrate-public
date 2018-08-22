clear all; close all; clc;

%% Select source data folder
folder = 'D:\Kieran\MATLAB_Scripts\LPA\LPA05\LPA05_rawData_IRIS\';
measurementNames = {'IRIS250.txt'; 'IRIS500.txt'; 'IRIS1000.txt'; 'IRIS2000.txt'; 'IRIS4000.txt'};
measurementValues = [250 500 1000 2000 4000]';
numWells = 8;

%% Set parameters
ampthresh = 0.7; % Percent of max intensity threshold for segmenting wells
sdthresh = 1.5; % Threshold for discarding unwanted data points from wells (ie, data captured when moving sensor to well)
minPeakDist = 0;

%% Import data
numMeasurements = length(measurementValues);
data = nan(100,numMeasurements);

for f = 1:numMeasurements
    rawData = readtable([folder measurementNames{f}]); rawData = wrev(rawData{:,2});    
    data(1:length(rawData),f) = rawData;
end

%% Measure light intensity per well
for i = 1:numMeasurements
    time = 1:length(data(:,i));
    subplot(numMeasurements,1,i); hold on;
    plot(time,data(:,i)); title([ regexprep(measurementNames{i},'.txt.*','') ' intensity ']); ylabel('Intensity (W)')
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
    if max(idx) ~= numWells
        disp(['Warning! More than ' num2str(numWells) ' wells detected for ' measurementNames{i}]);
    end
    
    %% Calculate average well intensities and exlude low outliers (captured when moving sensor to well)
    for j = 1:numWells
        m = nanmedian(wellIntensity(wellID==j));
        sd = nanstd(wellIntensity(wellID==j));
        wellIntensity(wellIntensity<(m-sdthresh*sd) & wellID==j) = nan;
        wellMean(i,j) = nanmean(wellIntensity(wellID==j));
        wellSD(i,j) = nanstd(wellIntensity(wellID==j));
    end
    
    plot(time,wellIntensity,'ro');
    
end

%% Calculate plate intensity per IRIS value
intensity = [zeros(8,1), wellMean'];
measurementValues = [0; measurementValues];
u = mean(intensity)';
sd = std(intensity)';

%% Fit data
[f, g] = fit(measurementValues, u,'poly1');
x = 0:1:max(measurementValues);
yfit = f.p1*x + f.p2;

%% Calculate IRIS needed for target light output
syms inputIRIS;
lightOutput= 100E-6;
eqn = f.p1*inputIRIS + f.p2;
inputIRIS = round(double(solve(eqn == lightOutput,inputIRIS)));
disp(['For target power = ' num2str(lightOutput*1E6) ' uW, inputIRIS = ' num2str(inputIRIS)])

%% Plot data
figure;
errorbar(measurementValues,u*1E6,sd*1E6,'ro','MarkerSize',2);
hold on;
plot(x,yfit*1E6);
title(['LPA response curve (R squared = ' num2str(g.rsquare) ')']); xlabel('IRIS values (AU)'); ylabel('Light intensity (uW)')
