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

%% Set calibration round and channel
%calibrationRound = 1;
UserAnswer = inputdlg('Calibration Round?')
calibrationRound = str2num(UserAnswer{1});
%8/23 Megan is lazy and wanted to just enter this data instead of changing
%the script

%% Select source data folders
%folder_ch1 = 'D:\Kieran\MATLAB_Scripts\LPA\LPA04\rawData_calibrate\Round4_channel1\';
%folder_ch2 = 'D:\Kieran\MATLAB_Scripts\LPA\LPA04\rawData_calibrate\Round4_channel2\';
%8/23 This now happens in formatData (MM)

%% Select output folder
%outputFolder = 'D:\Kieran\MATLAB_Scripts\LPA\LPA04\gcal';
outputFolder=uigetdir('','Select folder to store calibration files');
%8/23 Megan is lazy and wanted to just pick the output folder

%% Set parameters
ampthresh = 0.7; % Percent of max intensity threshold for segmenting wells
sdthresh = 1.5; % Threshold for discarding unwanted data points from wells (ie, data captured when moving sensor to well)
minPeakDist = 0;
maxCal = 255; % Initial max value of LED intensity

numRows = 4; % Rows A-D of 24 well plate
numColumns = 6; % Columns 1-6 of 24 well plate
numWells=numRows*numColumns;
rowNames = ['A'; 'B'; 'C'; 'D'];
channelsPerWell = 2;
%totalColumns = numColumns*channelsPerWell; %8/23 MM changed because new
%power meter software can read in all 24 wells
totalColumns=numWells*channelsPerWell;


%% Measure intensity from raw power meter data and assign to well
%folders = {folder_ch1, folder_ch2};
% 8/23 I changed this all because reading in data is happening in
% formatData function now so that only that file needs to be changed if the
% format of the data file (power meter) file changes MM

for c = 1:channelsPerWell
%     folder = folders{c}; 8/23 MM
    figure('Name', ['Channel ' num2str(c) ' well identification']);
    
    %% Import data from Thorlabs power meter
%8/23 MM Changed because all 24 wells can now be read using new power meter
%software.  So all the well data is now read into data (for one
%channel/LED) using the file formatData
%     A = readtable([folder 'A.txt']); A = wrev(A{:,2});
%     B = readtable([folder 'B.txt']); B = wrev(B{:,2});
%     C = readtable([folder 'C.txt']); C = wrev(C{:,2});
%     D = readtable([folder 'D.txt']); D = wrev(D{:,2});
    
%     data = nan(100,4);
%     data(1:length(A),1) = A;
%     data(1:length(B),2) = B;
%     data(1:length(C),3) = C;
%     data(1:length(D),4) = D;
      data=formatData();
    
    %% Measure light intensity per well
   % for i = 1:numRows
   i=1;   %8/23 MM Lazy fix due to no longer needing to lopp through rows. Should be changed when we're sure this is the way we'd like to implement 
  
   time = 1:length(data(:,1));
        subplot(4,1,1); hold on;
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
        if max(idx) ~= numColumns %8/23 NOTE: This warning probably needs to go bye-bye
            disp(['Warning! More than ' num2str(numColumns) ' wells detected in channel ' num2str(c) ' row ' rowNames(i)]);
        end
        
        %% Calculate average well intensities and exlude low outliers (captured when moving sensor to well)
        for j = 1:numWells
            m = nanmedian(wellIntensity(wellID==j));
            sd = nanstd(wellIntensity(wellID==j));
            wellIntensity(wellIntensity<(m-sdthresh*sd) & wellID==j) = nan;
            wellMean(c,i,j) = nanmean(wellIntensity(wellID==j));
            wellSD(c,i,j) = nanstd(wellIntensity(wellID==j));
        end
        
        plot(time,wellIntensity,'ro');
        
%end
end

%% Format measured intensity data for calibration
%rawIntensity = zeros(numRows, totalColumns);
rawIntensity=zeros(1,totalColumns);
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
%MM 8/23 Does this need to be output as 4X6 matrix to be read appropriately
%by the LPA software?--Doesn't seem to need to be in a matrix for
%appropriate LPA function.
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
    %heatmap(data{i});
    %8/23 MM changing to plot the linear data so it looks like the plate
        PlateFormat=vec2mat(data{i}(:),4);
        heatmap(PlateFormat)
    colorbar;
    title(titles{i});
end

%clearvars -except calibrationRound outputFolder maxCal numRows totalColumns wellMean wellSD