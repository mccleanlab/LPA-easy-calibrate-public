function Intensities=formatData()
% SUMMARY:
% This function combines the tracked results files in varargin (enter as
% many as you want, but more than 1).  Results files have the formed
% returned by trackCellsLocal or similar functions.
%
% INPUTS:
%  none
%
%
% OUTPUTS:
%  TBD
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
%% Open the appropriate data files

[FILENAME, PATHNAME, FILTERINDEX] = uigetfile('*.csv')
FID=fopen(strcat(PATHNAME,FILENAME));
Intensities=[];
while feof(FID)==0
    [C,position]= textscan(FID,'%f %s %s %f ','HeaderLines', 15);
    Intensities=[Intensities; C{4}];
end
fclose(FID);

end
