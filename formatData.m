function Intensities=formatData()
% SUMMARY:
% This function allows the user to choose (using uigetfile) an appropriate
% .csv file containing intensity measurements generation by the ThorLabs
% power meter. The intensity measurements in that file are returned in the
% variable Intensities.
%
% NOTE: This function is specific to the file format from the ThorLabs
% Optical Power Monitor v1.0.2149.55 software. If the format of the data
% file generated by the power meter changes this function may not work.
%
%
% INPUTS:
%  none
%
%
% OUTPUTS:
%  Intensities-an array containing the intensity measurements
%    
%
%    Written by Megan McClean, Ph.D.
%               University of Wisconsin-Madison
%               Department of Biomedical Engineering
%               1550 Engineering Drive ECB 3156
%               Madison, WI 53705
%               mmcclean@wisc.edu
%
%
%    Last revised on August 23, 2018
%
%% Open the appropriate data file
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile('*.csv');
FID=fopen(strcat(PATHNAME,FILENAME));
Intensities=[];

while feof(FID)==0 %keep scanning until the end of the FID file is reached
 [C,position]= textscan(FID,'%f %s %s %f');%,'HeaderLines', 15);
 if (feof(FID)==0) %check we aren't at the end of the file
        fseek(FID, position+1,'bof'); %advance to the next bit from where textscan stopped
 end
 Intensities=[Intensities; C{4}];
end
fclose(FID);
%Note: A loop is used here to avoid failure when the power meter file has
%an error, such as an intensity which is read as a random character, which
%happens every once in awhile. The error in format will cause textscan to
%abort, however, we just continue looping and scanning (and placing
%Intensities into the Intensities variable) until the end of the file is
%reached.

end
